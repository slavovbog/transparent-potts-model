using ExactPottsModel
Trajs = ExactPottsModel.Trajectories
using Random
using CairoMakie
using LaTeXStrings
CairoMakie.activate!()

empty!(Trajs.cache)
Trajs.read_records()

function plot_traj!(fig, q, n, T, J, n_steps, traj_id)
    config = (
        q = q, 
        n = n, 
        T = T, 
        J = J, 
        n_steps = n_steps, 
        burn_in = 0, 
        thin = 1,
        rng = Random.default_rng(),
    )
    outputs = Trajs.trajectory_outputs(config)
    @assert length(outputs) >= traj_id "traj_id is out of range"
    traj = outputs[traj_id].traj
    c = traj.c
    t = traj.t_sweeps

    ax = Axis(fig[1, 1]; xlabel = L"\text{MC time } t", ylabel = L"c_p(t)")
    ylims!(ax, -0.01, 1.01)

    for p in 1:q
        lines!(
            ax,
            t,
            c[p, :];
            ExactPottsModel.MakiePlotting.figure_style().traj_lines...,
            color = p,
            colorrange = (1, q),
        )
    end

    c_min = ExactPottsModel.global_argmin_1d_per_site_potential(
        q,
        Inf,
        T,
        J;
        tol = 1e-8,
        optimize_n_starts = 8,
    )
    length(c_min) > 1 && @warn("multiple minima found for q = $q, n = $n, T = $T, J = $J")
    c_min = first(c_min)
    hlines!(ax, [c_min], linestyle = :dash, color = :red)

    return fig
end

function figure6()
    q = 3
    n = 25
    T = 1
    Js = [(10, -5), (10, -15), (10, -40), (10, -21), (20, -40)]
    j1 = 1
    j2 = 3
    j3 = 5
    n_steps = 5000
    romannum(id) = begin
        id == 1 && return "I"
        id == 2 && return "II"
        id == 3 && return "III"
        id == 4 && return "IV"
        id == 5 && return "V"
        return "N/A"
    end

    fig = ExactPottsModel.MakiePlotting.generate_figure(size = (700, 250))

    fa = fig[1, 1]
    fb = fig[1, 2]
    fc = fig[1, 3]

    plot_traj!(fa, q, n, T, Js[j1], n_steps, 4)
    ax = content(fa[1, 1])
    ax.title = L"J_2=%$(Js[j1][1]),\ J_3=%$(Js[j1][2])\ \text{(%$(romannum(j1)))}"

    plot_traj!(fb, q, n, T, Js[j2], n_steps, 1)
    ax = content(fb[1, 1])
    ax.title = L"J_2=%$(Js[j2][1]),\ J_3=%$(Js[j2][2])\ \text{(%$(romannum(j2)))}"

    plot_traj!(fc, q, n, T, Js[j3], n_steps, 1)
    ax = content(fc[1, 1])
    ax.title = L"J_2=%$(Js[j3][1]),\ J_3=%$(Js[j3][2])\ \text{(%$(romannum(j3)))}"

    for (label, layout) in zip(["(a)", "(b)", "(c)"], [fa, fb, fc])
        Label(layout[1, 1, TopLeft()], label; padding = (0, 50, 5, 0))
    end

    return fig
end

function figure7()
    q = 3
    n = 25
    Js = [(10, -5), (10, -15), (10, -40), (10, -21), (20, -40)]
    j = 1
    T1 = 2
    T2 = 5
    n_steps = 5000
    romannum(id) = begin
        id == 1 && return "I"
        id == 2 && return "II"
        id == 3 && return "III"
        id == 4 && return "IV"
        id == 5 && return "V"
        return "N/A"
    end

    fig = ExactPottsModel.MakiePlotting.generate_figure(size = (700, 250))

    fa = fig[1, 1]
    fb = fig[1, 2]

    plot_traj!(fa, q, n, T1, Js[j], n_steps, 1)
    ax = content(fa[1, 1])
    ax.title = L"T=%$(T1)"
    # xlims!(ax, 0, 500)

    plot_traj!(fb, q, n, T2, Js[j], n_steps, 1)
    ax = content(fb[1, 1])
    ax.title = L"T=%$(T2)"
    # xlims!(ax, 0, 500)

    for (label, layout) in zip(["(a)", "(b)"], [fa, fb])
        Label(layout[1, 1, TopLeft()], label; padding = (0, 50, 5, 0))
    end

    Label(
        fig[0, :, Top()],
        L"J_2 = %$(Js[j][1]),\ J_3 = %$(Js[j][2])\ \text{(%$(romannum(j)))}";
        padding = (0, 0, 0, 0),
    )
    rowsize!(fig.layout, 0, Relative(0.01))

    return fig
end

FIGURES_DIR = joinpath(@__DIR__, "..", "figures")
mkpath(FIGURES_DIR)
save(joinpath(FIGURES_DIR, "figure6.pdf"), figure6());
save(joinpath(FIGURES_DIR, "figure7.pdf"), figure7());
