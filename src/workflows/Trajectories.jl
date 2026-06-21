# Batch generation and caching of MC concentration trajectories.
module Trajectories

using StableRNGs: StableRNG

using ..NumericalUtils
using ..Simulations
using ..ExactPottsModel

const cache = Dict()
const DATA_DIR = normpath(joinpath(@__DIR__, "..", "..", "data", "trajectories/"))

# Fixed salt so seeds are stable across Julia versions and platforms (via StableRNGs).
const TRAJECTORY_SEED_SALT = 42

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

function trajectory_seed(config, replicate_id::Int)
    @assert replicate_id >= 1 "replicate_id must be >= 1"
    h = hash((
        TRAJECTORY_SEED_SALT,
        config.q,
        config.n,
        config.T,
        config.J,
        config.n_steps,
        config.burn_in,
        config.thin,
        replicate_id,
    ))
    return Int(mod(h, Int128(2)^31 - 2)) + 1
end

trajectory_rng(seed::Integer) = StableRNG(seed)

function build_trajectory_input(config)
    @assert isfinite(config.n) "For trajectories n must be finite"
    replicate_id = get(config, :replicate_id, 1)
    seed = get(config, :seed, trajectory_seed(config, replicate_id))
    return (
        output_kind = :traj,
        model_params = model_params(config),
        n_steps = config.n_steps,
        burn_in = config.burn_in,
        thin = config.thin,
        replicate_id = replicate_id,
        seed = seed,
    )
end

function run_trajectory(input)
    model_params = input.model_params
    rng = trajectory_rng(input.seed)

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
        rng = rng,
    )

    return (traj = traj,)
end

function trajectory_output(config)
    input = build_trajectory_input(config)
    records = get_or_run!(input, run_trajectory)
    return only(records).output
end

function parallel_run(config, n_trajectories; batch_size)
    inputs = [
        build_trajectory_input((; config..., replicate_id = i)) for i in 1:n_trajectories
    ]

    cache_file = cache_path(inputs[1])
    @assert all(cache_path(input) == cache_file for input in inputs) "All inputs must share one cache file"

    ExactPottsModel.parallel_run!(cache_file, inputs, run_trajectory; batch_size = batch_size)
end

export trajectory_seed
export trajectory_output
export parallel_run

end # module
