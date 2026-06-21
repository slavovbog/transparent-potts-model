# Multistart optimization helpers and the implicit p*-model mapping.
module NumericalUtils

using Optim
using Random

@inline xlogx(x::Real) = x <= eps() ? 0.0 : float(x) * log(float(x))

function global_argmin_1d(f, domain::Tuple{<:Real,<:Real}; tol, optimize_n_starts)
    lower_bound, upper_bound = float(domain[1]), float(domain[end])
    start_points = collect(range(lower_bound + tol, upper_bound - tol, length = optimize_n_starts))
    candidates = Vector{Tuple{Float64,Float64}}()

    for x0 in start_points
        result = Optim.optimize(
            x -> f(first(x)),
            [lower_bound],
            [upper_bound],
            [x0],
        )

        trial_arg = result.minimizer[1]
        trial_value = result.minimum

        if all(abs(trial_arg - x) > tol for (x, _) in candidates)
            push!(candidates, (trial_arg, trial_value))
        end
    end
    isempty(candidates) && error("no minima found; check optimization parameters")

    global_min_value = minimum(val for (_, val) in candidates)
    global_min_args = [arg for (arg, val) in candidates if abs(val - global_min_value) <= tol]

    return global_min_args
end

function global_argmin_1d(f, domain; tol)
    f_values = f.(domain)
    global_min_value = minimum(f_values)
    global_min_args =
        [x for (x, fx) in zip(domain, f_values) if abs(fx - global_min_value) <= tol]
    return global_min_args
end

function global_argmin_nd(
    f,
    n_dims::Integer,
    domain::Tuple{<:Real,<:Real};
    tol,
    optimize_n_starts,
    rng::AbstractRNG = Random.default_rng(),
)
    lower_bound, upper_bound = float(domain[1]), float(domain[end])
    margin = max(tol, 100 * eps(Float64))
    span = upper_bound - lower_bound
    lb = fill(lower_bound, n_dims)
    ub = fill(upper_bound, n_dims)
    candidates = Tuple{Vector{Float64},Float64}[]

    for _ in 1:optimize_n_starts
        x0 = lb .+ margin .+ rand(rng, n_dims) .* (span - 2margin)
        result = Optim.optimize(f, lb, ub, x0, Optim.Fminbox(Optim.LBFGS()))
        trial_arg = collect(Float64, result.minimizer)
        trial_value = result.minimum
        if all(maximum(abs, trial_arg .- x) > tol for (x, _) in candidates)
            push!(candidates, (trial_arg, trial_value))
        end
    end
    isempty(candidates) && error("no minima found; check optimization parameters")

    global_min_value = minimum(val for (_, val) in candidates)
    global_min_args = [arg for (arg, val) in candidates if abs(val - global_min_value) <= tol]
    return global_min_args
end

function global_argmin_nd(f, n_dims::Integer, domain; tol)
    axes = ntuple(_ -> domain, n_dims)
    grid_points = vec(collect(Iterators.product(axes...)))
    f_values = [f(p) for p in grid_points]
    global_min_value = minimum(f_values)
    global_min_args =
        [p for (p, fv) in zip(grid_points, f_values) if abs(fv - global_min_value) <= tol]
    return global_min_args
end

function solve_implicit_Potts_pstar_equations(t_1, t_2, t_3; T, tol)
    # Map p*-coefficients (t_1, t_2, t_3) to Potts couplings at given q (Eqs. in the paper).
    J_3 = q -> 12t_3 * (q - 1)^2 / ((q - 1)^2 - 1)
    J_2 = q -> (q - 1) / q * (4t_2 - 12t_3 / ((q - 1)^2 - 1))

    squared_residual =
        q -> begin
            from_eq1 = (-T * log(q - 1) - J_2(q) / (q - 1) - J_3(q) / 2 / (q - 1)^2 - 2t_1)^2
            from_eq2 =
                (
                    -T * log(q - 1) - 4t_2 / q - (6t_3 * (1 - 2 / q)) / ((q - 1)^2 - 1) - 2t_1
                )^2
            @assert abs(from_eq1 - from_eq2) < tol "residuals are not equal: $from_eq1 vs $from_eq2"
            return from_eq1
        end

    # q > 2 and 1 < q < 2 can land in different intervals; search both.
    q_interval_above_2 = (2 + tol, 10)
    q_candidates_above_2 =
        global_argmin_1d(squared_residual, q_interval_above_2; tol = tol, optimize_n_starts = 10)
    q_candidates_above_2 = filter(q -> squared_residual(q) < tol, q_candidates_above_2)

    q_interval_between_1_and_2 = (1 + tol, 2 - tol)
    q_candidates_between_1_and_2 = global_argmin_1d(
        squared_residual,
        q_interval_between_1_and_2;
        tol = tol,
        optimize_n_starts = 10,
    )
    q_candidates_between_1_and_2 =
        filter(q -> squared_residual(q) < tol, q_candidates_between_1_and_2)

    q_candidates = vcat(q_candidates_above_2, q_candidates_between_1_and_2, [3])
    q = q_candidates[argmin(squared_residual.(q_candidates))]

    @info(
        "Found $(length(q_candidates)) solutions to the implicit Potts p* equations:",
        residuals = join(squared_residual.(q_candidates), ", "),
        q = join(q_candidates, ", "),
        J_2 = join(J_2.(q_candidates), ", "),
        J_3 = join(J_3.(q_candidates), ", "),
        q_minimal_residual = q,
        J_2_minimal_residual = J_2(q),
        J_3_minimal_residual = J_3(q),
        residual_minimal_residual = squared_residual(q),
    )

    return q, J_2(q), J_3(q)
end

export xlogx
export global_argmin_1d
export global_argmin_nd
export solve_implicit_Potts_pstar_equations

end # module
