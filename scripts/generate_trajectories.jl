using Distributed
@everywhere using ExactPottsModel
Traj = ExactPottsModel.Trajectories

for n in [10, 25],
    T in [0.8, 1, 1.5, 2, 5],
    J in [(10, -5), (10, -15), (10, -40), (10, -21), (20, -40)],
    n_steps in [5000],
    n_trajectories in [4]

    config = (
        q = 3,
        n = n,
        T = T,
        J = J,
        n_steps = n_steps,
        burn_in = 0,
        thin = 1
    )

    Traj.parallel_run(config, n_trajectories; batch_size = 1)
end
