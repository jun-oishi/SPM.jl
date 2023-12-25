module SPMPlots

using ..SPMCore

import Plots

"""
    # Parameters
    - surface: Surface
    - 他 : Plots.plotに準じるキーワード引数
    # Returns
    - Plots.Plotオブジェクト
"""
function heatmap(surface::SPMCore.Surface; kw...)::Plots.Plot
    x_mesh = collect(1:size(surface.matrix, 2)) .* surface.horz_resolution
    y_mesh = collect(1:size(surface.matrix, 1)) .* surface.horz_resolution
    return Plots.heatmap(x_mesh, y_mesh, surface.matrix, aspectratio=:equal; kw...)
end

"""
    # Summary
    - 表面形状を3次元でプロットする
    # Parameters
    - surface: Surface
    - 他 : Plots.plotに準じるキーワード引数
    # Returns
    - Plots.Plotオブジェクト
"""
function plot3dView(surface::SPMCore.Surface; kw...)::Plots.Plot
    x_mesh = collect(1:size(surface.matrix, 2)) .* surface.horz_resolution
    y_mesh = collect(1:size(surface.matrix, 1)) .* surface.horz_resolution
    return Plots.plot(x_mesh, y_mesh, surface.matrix, st=:surface; kw...)
end

"""
    # Summary
    - 直線上でのラインプロファイルを既存のプロットに追加する
    # Parameters
    - plot : Plots.Plotオブジェクト
    - surface: Surface
    - direction: :x or :y (水平方向または垂直方向)
    - i_slice: プロファイルを取得する列または行番号
    - 他 : Plots.plotに準じるキーワード引数
    # Returns
    - Plots.Plotオブジェクト
"""
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

"""
    # Summary
    - 直線上でのラインプロファイルを新規のプロットに描画する
    # Parameters
    - surface: Surface
    - direction: :x or :y (水平方向または垂直方向)
    - i_slice: プロファイルを取得する列または行番号
    - 他 : Plots.plotに準じるキーワード引数
    # Returns
    - Plots.Plotオブジェクト
"""
function plotProfile(
        surface::SPMCore.Surface, direction::Symbol, i_slice::Integer; kw...
)::Plots.Plot
    return plotProfile!(nothing, surface, direction, i_slice; kw...)
end

end # module Plot