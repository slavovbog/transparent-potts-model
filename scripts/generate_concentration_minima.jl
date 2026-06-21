using Distributed
@everywhere using ExactPottsModel
CMin = ExactPottsModel.ConcentrationMinima

# Phase-diagram grid; each (J_2, J_3) pair stores one c_min value in a single JLD2 file
data_J_2_range = -200:200
data_J_3_range = -200:200

for q in [1.1, 1.9, 2, 2.1, 3, 4, 10],
    n in [Inf],
    T in [1.0, 0.1, 10]

    CMin.parallel_run(
        (q = q, n = n, T = T, optimize_n_starts = 8, tol = 1e-8),
        data_J_2_range,
        data_J_3_range;
        batch_size = 50000,
    )
end
