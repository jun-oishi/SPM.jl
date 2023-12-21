module SPMPlots

using ..SPMCore

import Plots

"""
    # Parameters
    - surface: Surface
    - title: タイトル
    - legend: 凡例の表示, デフォルトはtrue
    # Returns
    - Plots.Plot
"""
function heatmap(surface::SPMCore.Surface; kw...)::Plots.Plot
    x_mesh = collect(1:size(surface.matrix, 2)) .* surface.horz_resolution
    y_mesh = collect(1:size(surface.matrix, 1)) .* surface.horz_resolution
    return Plots.heatmap(x_mesh, y_mesh, surface.matrix, aspectratio=:equal; kw...)
end

function plot3dView(surface::SPMCore.Surface; kw...)::Plots.Plot
    x_mesh = collect(1:size(surface.matrix, 2)) .* surface.horz_resolution
    y_mesh = collect(1:size(surface.matrix, 1)) .* surface.horz_resolution
    return Plots.plot(x_mesh, y_mesh, surface.matrix, st=:surface; kw...)
end

function plotProfile!(
        plot::Union{Plots.Plot,Nothing}, surface::SPMCore.Surface, direction::Symbol, i_slice::Integer; kw...
)::Plots.Plot
    if plot === nothing
        plot = Plots.plot()
    end
    if direction == :x
        x_mesh = collect(1:size(surface.matrix, 1)) .* surface.horz_resolution
        y_mesh = surface.matrix[:, i_slice]
    elseif direction == :y
        x_mesh = collect(1:size(surface.matrix, 2)) .* surface.horz_resolution
        y_mesh = surface.matrix[i_slice, :]
    else
        throw(ArgumentError("direction must be :x or :y"))
    end
    return Plots.plot!(plot, x_mesh, y_mesh; kw...)
end

function plotProfile(
        surface::SPMCore.Surface, direction::Symbol, i_slice::Integer; kw...
)::Plots.Plot
    return plotProfile!(nothing, surface, direction, i_slice; kw...)
end

end # module Plot