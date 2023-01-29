module Plot

using ..SPMCore
import Plots

function heatmap(
    surface::SPMCore.Surface;
    title::String="", legend=true
)::Plots.Plot
    x_mesh = collect(1:size(surface.data, 1)) .* surface.resolution * DEFAULT_UNIT
    y_mesh = collect(1:size(surface.data, 2)) .* surface.resolution * DEFAULT_UNIT
    ret = Plots.heatmap(
        x_mesh, y_mesh, surface.data .* DEFAULT_UNIT,
        title=title, aspect_ratio=:equal, legend=legend
    )
    return ret
end

function plotProfile(
        surface::SPMCore.Surface, direction::Symbol, i_slice::Integer;
        title::String="", label::String=""
)::Plots.Plot
    if direction == :x
        x_mesh = collect(1:size(surface.data, 1)) .* surface.resolution * DEFAULT_UNIT
        y_mesh = surface.data[:, i_slice] * DEFAULT_UNIT
    elseif direction == :y
        x_mesh = collect(1:size(surface.data, 2)) .* surface.resolution * DEFAULT_UNIT
        y_mesh = surface.data[i_slice, :] * DEFAULT_UNIT
    else
        throw(ArgumentError("direction must be :x or :y"))
    end
    return Plots.plot(x_mesh, y_mesh, title=title, label=label)
end

end