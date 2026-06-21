# Shared Makie theme, figure constructor, and plot-element defaults.
module MakiePlotting

using Makie

function generate_figure(; kwargs...)
    axis_theme = Theme(
        Axis = (
            xlabelpadding = 0,
            ylabelpadding = 0,
            xminorticks = IntervalsBetween(4),
            yminorticks = IntervalsBetween(4),
            xminorticksvisible = true,
            yminorticksvisible = true,
        ),
    )

    theme = merge(theme_latexfonts(), axis_theme)
    set_theme!(theme)

    return Figure(; figure_padding = 10, kwargs...)
end

function figure_style()
    return (
        c_min_heatmap = (
            colormap = reverse(cgrad(:RdYlBu)),
            nan_color = :grey,
        ),
        c_min_colorbar = (
            minorticksvisible = false,
        ),
        c_min_scatter = (
            color = :red,
            marker = :x,
            markersize = 10,
        ),
        c_min_scatter_finite_n = (
            marker = :star4,
            markersize = 15,
        ),
        xlogx_line = (
            color = :black,
            linewidth = 2,
            linestyle = :dash,
        ),
        traj_lines = (
            colormap = :plasma,
        ),
    )
end

end # module
