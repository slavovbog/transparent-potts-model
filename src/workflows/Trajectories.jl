# Batch generation and caching of MC concentration trajectories.
module Trajectories

using Random

using ..NumericalUtils
using ..Simulations
using ..ExactPottsModel

const cache = Dict()
const DATA_DIR = normpath(joinpath(@__DIR__, "..", "..", "data", "trajectories/"))
read_records(path::AbstractString) = ExactPottsModel.read_records!(cache, path)
read_records(input) = ExactPottsModel.read_record_file!(cache, cache_path(input))
read_records() = ExactPottsModel.read_record_dir!(cache, DATA_DIR)
write_records() = ExactPottsModel.write_records(cache)
get_records(input) = ExactPottsModel.get_records(cache, cache_path(input), input)
get_or_run!(input, run) = ExactPottsModel.get_or_run!(cache, cache_path(input), input, run)

function cache_path(input)
    stem = string(input.output_kind)
    n_steps = Int(input.n_steps)
    stem *= "_n_steps=$(n_steps)"
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

function build_trajectory_input(config)
    @assert isfinite(config.n) "For trajectories n must be finite"
    return (
        output_kind = :traj,
        model_params = model_params(config),
        n_steps = config.n_steps,
        burn_in = config.burn_in,
        thin = config.thin,
        rng = config.rng,
    )
end

function run_trajectory(input)
    model_params = input.model_params

    # All chains start fully polarized in color 1.
    traj = Simulations.generate_trajectory(
        model_params.q,
        model_params.n,
        model_params.T,
        model_params.J;
        spins_init = fill(1, model_params.n),
        n_steps = input.n_steps,
        burn_in = input.burn_in,
        thin = input.thin,
        rng = input.rng,
    )

    return (traj = traj,)
end

function trajectory_outputs(config)
    input = build_trajectory_input(config)
    records = get_or_run!(input, run_trajectory)
    return [record.output for record in records]
end

function parallel_run(config, n_trajectories; batch_size)
    inputs = [
        build_trajectory_input((; config..., rng = Random.default_rng())) for _ in 1:n_trajectories
    ]

    cache_file = cache_path(inputs[1])
    @assert all(cache_path(input) == cache_file for input in inputs) "All inputs must share one cache file"

    ExactPottsModel.parallel_run!(cache_file, inputs, run_trajectory; batch_size = batch_size)
end

end # module
