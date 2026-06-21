using ExactPottsModel
CMin = ExactPottsModel.ConcentrationMinima
using Random
using CairoMakie
using LaTeXStrings
CairoMakie.activate!()

empty!(CMin.cache)
CMin.read_records()

const data_J_2_range = -200:200
const data_J_3_range = -200:200

# Records in each cache file follow the same J_2 × J_3 grid ordering as parallel_run
function record_id(J_2, J_3)
    i = findfirst(data_J_2_range .== J_2)
    j = findfirst(data_J_3_range .== J_3)
    return LinearIndices((data_J_2_range, data_J_3_range))[i, j]
end

function read_c_min_cache(q, n, T, J)
    @assert length(J) == 2 "J must be a tuple of length 2"
    input = (
        output_kind = :c_min,
        model_params = (q = q, n = n, T = T),
    )
    path = CMin.cache_path(input)
    id = record_id(J[1], J[2])
    record = CMin.cache[path][id]

    @assert record.input.model_params.q == q &&
        record.input.model_params.n == n &&
        record.input.model_params.T == T &&
        record.input.model_params.J == J "Parameters do not match"
    return record.output.c_min
end

function plot_per_site_potential1d!(fig, q, n, T, Js::Vector)
    ax = Axis(
        fig[1, 1];
        xlabel = L"c",
        ylabel = L"f_{\infty,q=%$(q)}(c)",
        limits = (nothing, nothing, -0.75, 0.25),
        ygridvisible = false,
    )

    cs = 0:0.01:1
    for J in Js
        shift = ExactPottsModel.per_site_potential(0, q, n, T, J)
        f = c -> ExactPottsModel.per_site_potential(c, q, n, T, J) - shift
        lines!(ax, cs, f; label = L"J_2=%$(J[1]),\ J_3=%$(J[2])")

        c_min = read_c_min_cache(q, n, T, J)
        f_min = ExactPottsModel.per_site_potential.(c_min, q, n, T, Ref(J)) .- shift
        scatter!(ax, c_min, f_min; ExactPottsModel.MakiePlotting.figure_style().c_min_scatter...)
    end

    lines!(
        ax,
        cs,
        c -> ExactPottsModel.xlogx(c) + ExactPottsModel.xlogx(1 - c);
        label = L"c \log c + (1-c) \log (1-c)",
        ExactPottsModel.MakiePlotting.figure_style().xlogx_line...,
    )

    Legend(fig[1, 2], ax; framevisible = false)

    return fig
end

# NaN in heatmaps marks parameter points with multiple degenerate minima.
function J_2_J_3_heatmap!(fig, q, n, T, J_2_range, J_3_range)
    c_mins = fill(NaN, length(J_2_range), length(J_3_range))
    for i in eachindex(J_2_range), j in eachindex(J_3_range)
        c_min = read_c_min_cache(q, n, T, (J_2_range[i], J_3_range[j]))
        c_mins[i, j] = length(c_min) == 1 ? first(c_min) : NaN
    end

    ax = Axis(
        fig[1, 1];
        xlabel = L"J_2",
        ylabel = L"J_3",
        limits = (J_2_range[1], J_2_range[end], J_3_range[1], J_3_range[end]),
    )

    hm = heatmap!(
        ax,
        J_2_range,
        J_3_range,
        c_mins;
        colorrange = (0, 1),
        ExactPottsModel.MakiePlotting.figure_style().c_min_heatmap...,
    )

    ticks = [0, 1 / q, 1]
    labels = [L"0", L"1/q", L"1"]
    Colorbar(
        fig[1, 2],
        hm;
        ticks = (ticks, labels),
        label = L"c_{\min}",
        ExactPottsModel.MakiePlotting.figure_style().c_min_colorbar...,
    )

    return fig
end

function figure1()
    q = 2
    n = Inf
    T = 1
    Js = [(1, 5), (1, 2), (1, 0), (1, -2)]

    fig = ExactPottsModel.MakiePlotting.generate_figure(size = (800, 250))

    fa = fig[1, 1]
    fb = fig[1, 2]

    J_2_range = data_J_2_range
    J_3_range = data_J_3_range

    plot_per_site_potential1d!(fa, q, n, T, Js)

    J_2_J_3_heatmap!(fb, q, n, T, J_2_range, J_3_range)
    ax = content(fb[1, 1])
    lines!(
        ax,
        J_2_range,
        J_2 -> 4 * T - 2 * J_2;
        label = L"J_3 = 4T - 2J_2",
        color = :black,
    )
    lines!(
        ax,
        J_2_range,
        J_2 -> - 2 * J_2;
        label = L"J_3 = - 2J_2",
        ExactPottsModel.MakiePlotting.figure_style().xlogx_line...,
    )
    axislegend(ax)

    for (label, layout) in zip(["(a)", "(b)"], [fa, fb])
        Label(layout[1, 1, TopLeft()], label; padding = (0, 50, 5, 0))
    end

    colsize!(fig.layout, 2, Relative(0.35))

    return fig
end

function figure2()
    fig = ExactPottsModel.MakiePlotting.generate_figure(size = (250, 230))

    ax = Axis(
        fig[1, 1],
        xlabel = L"c",
        ylabel = L"f_{\infty}(c)",
        xgridvisible = false,
        ygridvisible = false,
    )
    cs = 0:0.01:1

    t = (-1.342, 1.871, -0.756)
    pstar_f =
        L -> begin
            -2 * sum(t[k] * L^k for k = 1:length(t)) +
            ExactPottsModel.xlogx(L) +
            ExactPottsModel.xlogx(1 - L)
        end

    n = Inf
    T = 1.0
    q, J_2, J_3 =
        ExactPottsModel.solve_implicit_Potts_pstar_equations(t...; T = T, tol = 1e-5)
    J = (J_2, J_3)
    shift = ExactPottsModel.per_site_potential(0, q, n, T, J) - pstar_f(0)
    potts_f = c -> ExactPottsModel.per_site_potential(c, q, n, T, J) - shift

    lines!(ax, cs, potts_f; color = :black, linewidth = 3, label = "Potts model")
    lines!(
        ax,
        cs,
        pstar_f;
        color = :red,
        linestyle = :dash,
        label = latexstring("\$p\$-star model"),
    )

    axislegend(ax, position = :lt)

    return fig
end

function plot_per_site_potential3d!(fig, n, T, J)
    q = 3

    ax = Axis3(
        fig[1, 1],
        xlabel = L"c_1",
        ylabel = L"c_2",
        zlabel = L"F_\infty",
        title = length(J) == 2 ? L"J_2=%$(J[1]),\ J_3=%$(J[2])" : "J=%$(J)",
        aspect = :equal,
        elevation = 0.1π,
        azimuth = 1.35π,
        xlabeloffset = 30,
        ylabeloffset = 30,
        zlabeloffset = 50,
    )
    zlims!(ax, -3.7, 1.5)

    shift = ExactPottsModel.per_site_potential((1, 0, 0), n, T, J)
    
    n_points = 300
    x_range = range(0, 1, n_points)
    y_range = range(0, 1, n_points)
    x_grid = [x for x in x_range, _ in y_range]
    y_grid = [y for _ in x_range, y in y_range]
    z_grid = [
        ExactPottsModel.per_site_potential((x, y, 1 - x - y), n, T, J) - shift for
        x in x_range, y in y_range
    ]
    valid_mask = x_grid .+ y_grid .< 0.99
    x_valid = x_grid[valid_mask]
    y_valid = y_grid[valid_mask]
    z_valid = z_grid[valid_mask]
    surface!(
        ax,
        x_valid,
        y_valid,
        z_valid,
        colormap = cgrad(:deep),
        shading = NoShading,
        rasterize = 5,
    )

    contour!(
        ax,
        x_valid,
        y_valid,
        z_valid,
        levels = 20,
        colormap = cgrad(:deep),
        transformation = (:xy, -3.5),
        # colorrange = (-3.5, 1.5),
    )
    # zlims!(ax, minimum(z_valid) - 0.01, maximum(z_valid) + 0.01)

    c_min = first(read_c_min_cache(q, n, T, J))

    conf_1 = (c_min, (1 - c_min) / 2, (1 - c_min) / 2)
    conf_2 = ((1 - c_min) / 2, c_min, (1 - c_min) / 2)
    conf_3 = ((1 - c_min) / 2, (1 - c_min) / 2, c_min)
    xs = [cs[1] for cs in (conf_1, conf_2, conf_3)]
    ys = [cs[2] for cs in (conf_1, conf_2, conf_3)]
    zs = [
        ExactPottsModel.per_site_potential(cs, n, T, J) - shift for
        cs in (conf_1, conf_2, conf_3)
    ]
    scatter!(ax, xs, ys, zs; ExactPottsModel.MakiePlotting.figure_style().c_min_scatter...)

    return fig
end

function figure3()
    q = 3
    n = Inf
    T = 1

    fig =
        ExactPottsModel.MakiePlotting.generate_figure(size = (900, 500), figure_padding = (20, 0, 10, 10))

    f1 = fig[2, 1] = GridLayout()
    f2 = fig[2, 2] = GridLayout()
    f3 = fig[2, 3] = GridLayout()
    f4 = fig[2, 4] = GridLayout()
    f5 = fig[2, 5] = GridLayout()
    row1 = fig[1, :] = GridLayout()
    f0 = row1[1, 1] = GridLayout()

    J_2_J_3_heatmap!(f0, q, n, T, 0:25, -50:0)
    plot_per_site_potential3d!(f1, n, T, (10, -5))
    plot_per_site_potential3d!(f2, n, T, (10, -15))
    plot_per_site_potential3d!(f3, n, T, (10, -40))
    plot_per_site_potential3d!(f4, n, T, (10, -21))
    plot_per_site_potential3d!(f5, n, T, (20, -40))

    rowsize!(fig.layout, 1, Relative(0.7))
    colsize!(row1, 1, Aspect(1, 1.2))

    return fig
end

function figure4()
    qs = [1.1, 1.9, 2, 2.1, 4, 10]
    n = Inf
    T = 1

    fig = ExactPottsModel.MakiePlotting.generate_figure(size = (900, 450))

    for id in eachindex(qs)
        q = qs[id]
        i = (id - 1) ÷ 3 + 1
        j = (id - 1) % 3 + 1

        fij = fig[i, j] = GridLayout()
        J_2_J_3_heatmap!(fij, q, n, T, data_J_2_range, data_J_3_range)
        ax = content(fij[1, 1])
        ax.title = L"q=%$(q)"
    end

    return fig
end

figure4()

function plot_per_site_potential1d_with_inf_n!(fig, q, finite_n, T, Js, show_flags;)
    ax = Axis(
        fig[1, 1];
        xlabel = L"c",
        ylabel = L"f_n(c)",
        limits = (nothing, nothing, nothing, 2.5),
        ygridvisible = false,
    )

    romannum(id) = begin
        id == 1 && return "I"
        id == 2 && return "II"
        id == 3 && return "III"
        id == 4 && return "IV"
        id == 5 && return "V"
        return "N/A"
    end
    for (id, J, show_flag) in zip(eachindex(Js), Js, show_flags)
        if show_flag
            shift = ExactPottsModel.per_site_potential(0, q, finite_n, T, J)
            f = c -> ExactPottsModel.per_site_potential(c, q, finite_n, T, J) - shift
            scatterlines!(
                ax,
                0:(1/finite_n):1,
                f;
                color = Cycled(id),
            )

            shift_inf = ExactPottsModel.per_site_potential(0, q, Inf, T, J)
            f_inf = c -> ExactPottsModel.per_site_potential(c, q, Inf, T, J) - shift_inf
            lines!(
                ax,
                0:0.01:1,
                f_inf;
                color = Cycled(id),
                linestyle = :dash,
            )

            band!(
                ax,
                0:0.01:1,
                f_inf.(0:0.01:1),
                f.(0:0.01:1);
                color = Cycled(id),
                alpha = 0.2,
            )

            c_min = read_c_min_cache(q, finite_n, T, J)
            f_min =
                ExactPottsModel.per_site_potential.(c_min, q, finite_n, T, Ref(J)) .- shift
            scatter!(
                ax,
                c_min,
                f_min;
                ExactPottsModel.MakiePlotting.figure_style().c_min_scatter...,
                color = Cycled(id),
                marker = :star4,
                markersize = 15,
            )

            c_min_inf = read_c_min_cache(q, Inf, T, J)
            f_min_inf =
                ExactPottsModel.per_site_potential.(c_min_inf, q, Inf, T, Ref(J)) .-
                shift_inf
            scatter!(
                ax,
                c_min_inf,
                f_min_inf;
                ExactPottsModel.MakiePlotting.figure_style().c_min_scatter...,
            )
        end
    end

    legend_elems = []
    legend_labels = []
    for (id, J, show_flag) in zip(eachindex(Js), Js, show_flags)
        if show_flag
            push!(
                legend_elems,
                PolyElement(color = Cycled(id),)
            )
            push!(
                legend_labels,
                L"J_2=%$(J[1]),\ J_3=%$(J[2])\ \text{(%$(romannum(id)))}"
            )
        end
    end
    axislegend(
        ax,
        legend_elems,
        legend_labels;
        position = :lt,
    )

    return fig
end

function plot_c_min_vs_n!(fig, q, ns, T, Js)
    ax = Axis(
        fig[1, 1],
        xlabel = L"n",
        ylabel = L"c_{\min}",
        limits = (ns[1], ns[end], -0.1, 1.1),
        yminorticksvisible = false,
    )

    romannum(id) = begin
        id == 1 && return "I"
        id == 2 && return "II"
        id == 3 && return "III"
        id == 4 && return "IV"
        id == 5 && return "V"
        return "N/A"
    end
    ax_tick_vals = []
    ax_tick_labels = []
    legend_elems = []
    legend_labels = []
    for (id, J) in enumerate(Js)
        c_mins = []
        for n in ns
            c_min = ExactPottsModel.global_argmin_1d_per_site_potential(q, n, T, J, tol = 1e-8)
            length(c_min) > 1 &&
                @warn("multiple minima found for q = $q, n = $n, T = $T, J = $J")
            push!(c_mins, first(c_min))
        end
        scatterlines!(
            ax,
            ns,
            c_mins;
            ExactPottsModel.MakiePlotting.figure_style().c_min_scatter_finite_n...,
        )
        push!(legend_elems, LineElement(color = Cycled(id), linewidth = 3))
        push!(legend_labels, L"J_2=%$(J[1]),\ J_3=%$(J[2])\ \text{(%$(romannum(id)))}")

        c_min_inf = read_c_min_cache(q, Inf, T, J)
        length(c_min_inf) > 1 &&
            @warn("multiple minima found for q = $q, n = Inf, T = $T, J = $J")
        c_min_inf = first(c_min_inf)
        hlines!(ax, [c_min_inf], linestyle = :dash)
        push!(ax_tick_vals, round(c_min_inf, digits = 2))
        push!(ax_tick_labels, begin
            if c_min_inf ≈ 1 / q
                L"1/q"
            elseif c_min_inf ≈ 0.0
                L"0"
            elseif c_min_inf ≈ 1.0
                L"1"
            else
                L"%$(round(c_min_inf, digits = 2))"
            end
        end)
    end
    ax.yticks = (ax_tick_vals, ax_tick_labels)

    Legend(fig[1, 2], legend_elems, legend_labels; framevisible = false)

    return fig
end

function figure5()
    q = 3
    n = 10
    T = 1
    Js = [(10, -5), (10, -15), (10, -40), (10, -21), (20, -40)]
    show_flags = [true, false, true, false, true]

    fig = ExactPottsModel.MakiePlotting.generate_figure(size = (800, 300))

    fa = fig[1, 1] = GridLayout()
    fb = fig[1, 2] = GridLayout()

    plot_per_site_potential1d_with_inf_n!(fa, q, n, T, Js, show_flags)

    plot_c_min_vs_n!(fb, q, 2:50, T, Js)

    for (label, layout) in zip(["(a)", "(b)"], [fa, fb])
        Label(layout[1, 1, TopLeft()], label; padding = (0, 30, 5, 0))
    end

    colsize!(fig.layout, 1, Relative(0.35))
    colsize!(fb, 1, Relative(0.59))

    return fig
end

FIGURES_DIR = joinpath(@__DIR__, "..", "figures")
mkpath(FIGURES_DIR)
save(joinpath(FIGURES_DIR, "figure1.pdf"), figure1());
save(joinpath(FIGURES_DIR, "figure2.pdf"), figure2());
save(joinpath(FIGURES_DIR, "figure3_0.pdf"), figure3());
save(joinpath(FIGURES_DIR, "figure4.pdf"), figure4());
save(joinpath(FIGURES_DIR, "figure5.pdf"), figure5());
