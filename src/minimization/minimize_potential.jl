# Global minima of the symmetric per-site potential f_n(c) (1D) and f_n(c) (full qD).
using ..NumericalUtils
using ..ExactPottsFormulas

function global_argmin_1d_per_site_potential(
    q, n, T, J;
    check = true,
    tol,
    optimize_n_starts = nothing,
)
    potential(c) = per_site_potential(c, q, n, T, J; check = check)

    if isinf(n)
        # Thermodynamic limit: c is continuous on [0, 1]; need multistart search.
        @assert !isnothing(optimize_n_starts) "optimize_n_starts must be set when n = Inf"
        return global_argmin_1d(potential, (0, 1); tol = tol, optimize_n_starts = optimize_n_starts)
    else
        # Finite n: c = 0, 1/n, …, 1 is an exact grid.
        c_grid = 0:(1 / n):1
        return global_argmin_1d(potential, c_grid; tol = tol)
    end
end

function global_argmin_per_site_potential(
    q, n, T, J;
    tol,
    optimize_n_starts = nothing,
    rng = nothing,
)
    @assert isinteger(q) "q must be integer for the full potential"
    potential(cs) = per_site_potential(cs, n, T, J)

    if isinf(n)
        @assert !isnothing(optimize_n_starts) "optimize_n_starts must be set when n = Inf"
        @assert !isnothing(rng) "rng must be set when n = Inf"
        return global_argmin_nd(
            potential,
            q,
            (0, 1);
            tol = tol,
            optimize_n_starts = optimize_n_starts,
            rng = rng,
        )
    else
        c_grid = 0:(1 / n):1
        return global_argmin_nd(potential, q, c_grid; tol = tol)
    end
end
