# Hamiltonian and per-site potential f_n (notation as in the manuscript: q, n, r, J_k, n_p, c_p).
module ExactPottsFormulas

using ..NumericalUtils

using SpecialFunctions: gamma, loggamma
using Combinatorics: combinations, partitions

function Hamiltonian(spins, J, h = nothing)
    n = length(spins)
    r = length(J) + 1

    H = 0.0
    for k in 2:r
        J_k = J[k - 1]

        # Count k-tuples of sites with the same spin.
        n_same_spin_k_tuples = 0.0
        for site_indices in combinations(1:n, k)
            s_ref = spins[first(site_indices)]
            all_same = all(spins[i] == s_ref for i in site_indices)
            n_same_spin_k_tuples += all_same ? 1 : 0
        end

        H -= J_k / n^(k - 1) * n_same_spin_k_tuples
    end

    if !isnothing(h)
        q = length(h)
        @assert all(s ∈ 1:q for s in spins) "each spin must be in 1:q"

        H -= sum(h[s] for s in spins)
    end

    return H
end

function Hamiltonian(rep::Symbol, ns, J, h = nothing)
    if rep == :occupation_numbers
        q = length(ns)
        n = sum(ns)
        r = length(J) + 1

        H = 0.0
        for k in 2:r
            J_k = J[k - 1]
            # sum_p binom(n_p, k): all same-color k-tuples on the complete graph.
            n_k_tuples = sum(binomial(n_p, k) for n_p in ns)
            H -= J_k / n^(k - 1) * n_k_tuples
        end

        if !isnothing(h)
            @assert length(h) == q "length(h) must equal q"

            H -= sum(h[p] * ns[p] for p in 1:q)
        end

        return H
    else
        error("$rep is not a valid representation")
    end
end

function per_site_potential(ns, T, J)
    n = sum(ns)

    H = Hamiltonian(:occupation_numbers, ns, J)
    # S = log(n!) - sum_p log(n_p!).
    S = loggamma(1 + n) - sum(loggamma(1 + n_p) for n_p in ns)

    return (H - T * S) / n
end

function per_site_potential(cs, n, T, J)
    @assert sum(cs) ≈ 1.0 "concentrations c_p must sum to 1"
    isfinite(n) && return per_site_potential(n .* cs, T, J)

    q = length(cs)
    r = length(J) + 1

    H = 0.0
    for p in 1:q
        c_p = cs[p]
        for k in 2:r
            J_k = J[k - 1]
            H -= J_k / gamma(1 + k) * c_p^k
        end
    end

    S = -sum(xlogx.(cs))

    return H - T * S
end

function per_site_potential(c, q, n, T, J; check = true)
    !isinf(n) && !isinteger(q) && error("if n is finite, q must be an integer")

    if isinf(n)
        c = float(c)
        r = length(J) + 1

        # Symmetric ansatz: c_1 = c, c_p = (1-c)/(q-1) for p > 1.
        U_part = 0.0
        for k in 2:r
            J_k = J[k - 1]
            U_part -= J_k / gamma(1 + k) * ((1 - c)^k / (q - 1)^(k - 1) + c^k)
        end
        U_part += c * T * log(q - 1)
        U_part -= T * log(q - 1)

        S = -xlogx(c) - xlogx(1 - c)

        f_n = U_part - T * S

        if check && isinteger(q) && 1e-8 < c < 1 - 1e-8
            q = Int(q)
            cs = [c; fill((1 - c) / (q - 1), q - 1)]
            f_n_from_cs = per_site_potential(cs, n, T, J)

            @assert abs(f_n - f_n_from_cs) < 1e-8 "symmetric and vector forms disagree for c = $c, q = $q, n = $n, T = $T, J = $J"
        end
    else
        q = Int(q)
        cs = [c; fill((1 - c) / (q - 1), q - 1)]
        f_n = per_site_potential(cs, n, T, J)
    end

    return f_n
end

function check_Hamiltonian_representations(q, n, J, h = nothing; tol::Real)
    n = Int(n)
    q = Int(q)

    spin_ranges = ntuple(_ -> 1:q, n)

    max_abs_diff = 0.0
    n_checked = 0

    for spin_config in Iterators.product(spin_ranges...)
        spins = collect(spin_config)

        H_spins = Hamiltonian(spins, J, h)

        ns = zeros(Int, q)
        for s in spins
            ns[s] += 1
        end

        H_ns = Hamiltonian(:occupation_numbers, ns, J, h)

        diff = abs(H_spins - H_ns)
        max_abs_diff = max(max_abs_diff, diff)
        n_checked += 1

        if diff > tol
            error(
                "Hamiltonian representations disagree (abs diff = $diff, tol = $tol) for configuration $spins",
            )
        end
    end

    @info(
        "check_Hamiltonian_representations passed",
        n = n,
        q = q,
        J = J,
        tol = tol,
        n_checked = n_checked,
        max_abs_diff = max_abs_diff,
    )

    return true
end

function check_Hamiltonian_representations(; tol::Real)
    for q in [2, 3, 4], n in [2, 3, 4, 5], r in 2:4
        J = randn(r - 1)
        check_Hamiltonian_representations(q, n, J; tol = tol)
    end
end

export Hamiltonian
export per_site_potential
export check_Hamiltonian_representations

end # module
