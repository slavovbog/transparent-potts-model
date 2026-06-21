# Metropolis MC in the occupation-number representation (single-spin flips).
module Simulations

using Combinatorics: binomial
using Random

using ..ExactPottsFormulas: Hamiltonian

function flip_energy_change(ns, a, b, J, n)
    r = length(J) + 1
    n_a = ns[a]
    n_b = ns[b]

    ΔH = 0.0
    for k in 2:r
        J_k = J[k - 1]

        # Change in sum_p binom(n_p, k) when one spin moves a -> b.
        Δ_n_k_tuples =
            (binomial(n_a - 1, k) - binomial(n_a, k)) +
            (binomial(n_b + 1, k) - binomial(n_b, k))
        ΔH -= J_k / n^(k - 1) * Δ_n_k_tuples
    end

    return ΔH
end

function generate_trajectory(q, n, T, J; spins_init, n_steps, burn_in, thin, rng)
    @assert isinteger(q) "q must be an integer"
    @assert isfinite(n) "n must be finite"
    @assert isinteger(n) "n must be an integer"
    @assert T >= 0 "T must be nonnegative"
    @assert q >= 2 "q must be >= 2"
    @assert n >= 1 "n must be positive"
    @assert n_steps >= 1 "n_steps must be >= 1"
    @assert burn_in >= 0 "burn_in must be >= 0"
    @assert thin >= 1 "thin must be >= 1"
    @assert length(spins_init) == n "spins_init must have length n"
    @assert all(s ∈ 1:q for s in spins_init) "spins_init must be in 1:q"

    ns = zeros(Int, q)
    spins = copy(spins_init)
    for s in spins
        ns[s] += 1
    end

    H = Hamiltonian(:occupation_numbers, ns, J)

    n_recorded = n_steps <= burn_in ? 0 : (n_steps - burn_in) ÷ thin
    c = Matrix{Float64}(undef, q, n_recorded)
    t_sweeps = Vector{Int}(undef, n_recorded)

    accept_count = 0
    n_attempts = 0
    record_idx = 0

    for sweep in 1:n_steps
        for _ in 1:n
            i = rand(rng, 1:n)
            a = spins[i]

            # Uniform proposal among the other q-1 colors.
            t = rand(rng, 1:(q - 1))
            b = t < a ? t : t + 1

            ΔH = flip_energy_change(ns, a, b, J, n)

            accept = false
            if T == 0
                accept = ΔH <= 0
            elseif ΔH <= 0
                accept = true
            else
                accept = log(rand(rng)) < (-ΔH / T)
            end

            n_attempts += 1
            if accept
                spins[i] = b
                ns[a] -= 1
                ns[b] += 1
                H += ΔH
                accept_count += 1
            end
        end

        if sweep > burn_in && ((sweep - burn_in) % thin == 0)
            record_idx += 1
            c[:, record_idx] .= Float64.(ns) ./ n
            t_sweeps[record_idx] = sweep
        end
    end

    accept_rate = accept_count / max(1, n_attempts)
    return (c = c, t_sweeps = t_sweeps, accept_rate = accept_rate)
end

export generate_trajectory

end # module
