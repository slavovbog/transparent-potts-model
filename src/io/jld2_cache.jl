# Read / write JLD2 record files

function read_records(path)
    if !isfile(path)
        @warn("File $path does not exist")
        return []
    end

    records = []
    jldopen(path, "r") do file
        # Keys are "1", "2", …; sort numerically to preserve grid order.
        sorted_keys = sort(collect(keys(file)), by = k -> parse(Int, k))
        for key in sorted_keys
            push!(records, file[key])
        end
    end

    @info "Loaded $(length(records)) entries from $path"
    return records
end

function read_record_file!(cache, path)
    path ∈ keys(cache) && @warn("Read data will be appended")
    path ∉ keys(cache) && (cache[path] = [])
    append!(cache[path], read_records(path))
end

function read_record_dir!(cache, dir)
    for filename in readdir(dir)
        endswith(filename, ".jld2") || continue
        read_records!(cache, normpath(joinpath(dir, filename)))
    end
end

function read_records!(cache, path)
    isdir(path) && return read_record_dir!(cache, path)
    isfile(path) && return read_record_file!(cache, path)
end

function write_records!(path, records)
    @assert !isfile(path) "File $path already exists"

    record_index = 1
    jldopen(path, "w") do file
        for record in records
            file["$record_index"] = record
            record_index += 1
        end
    end

    @info "Wrote $(record_index - 1) data points to $path"
end

function write_records(cache)
    for (path, records) in cache
        write_records!(path, records)
    end
end

function cache_record!(cache, path, record)
    push!(get!(cache, path, []), record)
end

function get_records(cache, path, input)
    records = get(cache, path, [])

    matching_records = [record for record in records if record.input == input]

    length(matching_records) == 0 && @info("No cached entry for input $input in $path")
    length(matching_records) > 1 &&
        @warn("Found $(length(matching_records)) cached entries for input $input in $path")

    return matching_records
end

function get_or_run!(cache, path, input, run)
    matching_records = get_records(cache, path, input)
    if length(matching_records) == 0
        record = (input = input, output = run(input))
        cache_record!(cache, path, record)
        return [record]
    else
        return matching_records
    end
end

# Distributed batch jobs with incremental JLD2 append

function parallel_run!(path, inputs, run; batch_size::Integer)
    @assert batch_size ≥ 1 "batch_size must be at least 1"

    next_record_index = 0
    if isfile(path)
        @warn "Cache file $path already exists. Will append to it."
        jldopen(path, "r") do file
            for key in keys(file)
                record_index = parse(Int, key)
                record = file[key]
                @assert :input ∈ keys(record) && :output ∈ keys(record) "At least one record is missing an input or output"
                next_record_index = max(next_record_index, record_index)
            end
        end
    end
    next_record_index += 1

    n_total = length(inputs)

    @info "Parallel run: $n_total inputs ($(nworkers()) workers)"

    t_start = time()
    n_completed = 0
    n_failed = 0

    for batch_start = 1:batch_size:n_total
        batch_end = min(batch_start + batch_size - 1, n_total)
        batch_inputs = @view inputs[batch_start:batch_end]

        t_batch = time()
        outputs = pmap(run, batch_inputs; on_error = identity)
        dt_batch = time() - t_batch

        jldopen(path, isfile(path) ? "a" : "w") do file
            for (input, output) in zip(batch_inputs, outputs)
                if output isa Exception
                    @warn "Run failed for an input" exception = output
                    n_failed += 1
                    continue
                end
                file["$next_record_index"] = (input = input, output = output)
                next_record_index += 1
                n_completed += 1
            end
        end

        elapsed = time() - t_start
        rate = n_completed / elapsed
        remaining = (n_total - batch_end) / max(rate, eps())
        @info "Progress: $n_completed/$n_total " *
              "(batch in $(round(dt_batch; digits=1))s, " *
              "ETA $(round(remaining / 60; digits=1)) min)"
        flush(stdout)
    end

    elapsed_total = time() - t_start
    @info "Parallel run complete: $n_completed succeeded, $n_failed failed " *
          "in $(round(elapsed_total / 60; digits=1)) min"
end
