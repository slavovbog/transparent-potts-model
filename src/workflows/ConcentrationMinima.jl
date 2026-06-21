# Batch computation and caching of c_min(J_2, J_3) on the phase-diagram grid.
module ConcentrationMinima

using ..NumericalUtils
using ..ExactPottsModel

const cache = Dict()
const DATA_DIR = normpath(joinpath(@__DIR__, "..", "..", "data", "c_min_J_2_J_3/"))
read_records(path::AbstractString) = ExactPottsModel.read_records!(cache, path)
read_records(input) = ExactPottsModel.read_record_file!(cache, cache_path(input))
read_records() = ExactPottsModel.read_record_dir!(cache, DATA_DIR)
write_records() = ExactPottsModel.write_records(cache)
get_records(input) = ExactPottsModel.get_records(cache, cache_path(input), input)
get_or_run!(input, run) = ExactPottsModel.get_or_run!(cache, cache_path(input), input, run)

function cache_path(input)
    stem = string(input.output_kind)
    q_val = Float64(input.model_params.q)
    n_val = Float64(input.model_params.n)
    T_val = Float64(input.model_params.T)
    stem *= "_q=$(q_val)_n=$(n_val)_T=$(T_val)"
    return normpath(joinpath(DATA_DIR, "$(stem).jld2"))
end

function model_params(config)
    return (
        q = config.q,
        n = config.n,
        T = config.T,
        J = config.J,
    )
end

function build_c_min_input(config)
    if isinf(config.n)
        return (
            output_kind = :c_min,
            model_params = model_params(config),
            tol = config.tol,
            optimize_n_starts = config.optimize_n_starts,
        )
    else
        return (
            output_kind = :c_min,
            model_params = model_params(config),
            tol = config.tol,
        )
    end
end

function compute_c_min(input)
    model_params = input.model_params

    if isinf(model_params.n)
        c_min = ExactPottsModel.global_argmin_1d_per_site_potential(
            model_params.q,
            model_params.n,
            model_params.T,
            model_params.J;
            check = false,
            tol = input.tol,
            optimize_n_starts = input.optimize_n_starts,
        )
    else
        c_min = ExactPottsModel.global_argmin_1d_per_site_potential(
            model_params.q,
            model_params.n,
            model_params.T,
            model_params.J;
            check = false,
            tol = input.tol,
        )
    end

    return (c_min = c_min,)
end

function c_min_output(config)
    input = build_c_min_input(config)
    records = get_or_run!(input, compute_c_min)
    return only(records).output
end

function parallel_run(config, J_2_values, J_3_values; batch_size)
    inputs = [
        build_c_min_input((; config..., J = (J_2, J_3))) for
        J_2 in J_2_values, J_3 in J_3_values
    ]

    cache_file = cache_path(inputs[1])
    @assert all(cache_path(input) == cache_file for input in inputs) "All inputs must share one cache file"

    ExactPottsModel.parallel_run!(cache_file, inputs, compute_c_min; batch_size = batch_size)
end

end # module
