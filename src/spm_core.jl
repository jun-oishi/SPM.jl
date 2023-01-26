"""
STMデータを扱うためのモジュール
"""
module SPMCore

import DelimitedFiles.readdlm
using Plots

export Plots
export Surface, show
export Image, loadImage, extract, filter, horz_diff, vert_diff
export Tip

########################################################
### utilities ##########################################
########################################################

is_square(x::AbstractMatrix) = size(x, 1) == size(x, 2)

throw_if(x::Bool, msg::String) = x && throw(ArgumentError(msg))

abstract type Surface end
(s::Surface)(x::Integer, y::Integer) = s.data[y, x]
(s::Surface)(xrange::UnitRange, yrange::UnitRange) = s.data[yrange, xrange]

function showHeatmap(
    surface::Surface;
    title::String="", legend=false
)::Plots.Plot
    pyplot()
    x_mesh = collect(1:size(surface.data, 1)) .* surface.resolution / 10.0
    y_mesh = collect(1:size(surface.data, 2)) .* surface.resolution / 10.0
    ret = heatmap(
        x_mesh, y_mesh, surface.data, xlabel="x (nm)", ylabel="y (nm)",
        title=title, aspect_ratio=:equal, legend=legend
    )
    return ret
end

function addTompgraph!(
        plot::Plots.Plot, surface::Surface, direction::Symbol, i_slice::Integer;
        title::String="", label::String=""
)::Plots.Plot
    pyplot()
    if direction == :x
        x_mesh = collect(1:size(surface.data, 1)) .* surface.resolution / 10.0
        y_mesh = surface.data[:, i_slice]
    elseif direction == :y
        x_mesh = collect(1:size(surface.data, 2)) .* surface.resolution / 10.0
        y_mesh = surface.data[i_slice, :]
    else
        throw(ArgumentError("direction must be :x or :y"))
    end
    ret = plot!(
        plot, x_mesh, y_mesh, xlabel="horz. dist. (nm)", ylabel="height (nm)",
        title=title, label=label
    )
    return ret
end

########################################################
### Image ##############################################
########################################################

"""
    # Summary
    - 画像データを扱うための型
    # Fields
    - data: 画像データ
    - resolution: 解像度 (AA / px)
    - width: 幅 (AA)
    - height: 高さ (AA)
    # Constructer Parameters

"""
struct Image <: Surface
    data::Matrix
    resolution::AbstractFloat  # \AA / px
    width::AbstractFloat       # \AA
    height::AbstractFloat      # \AA

    """
        # Parameters
        - data: 画像データ
        - resolution: 解像度 (AA / px)
    """
    function Image(data::Matrix; resolution::AbstractFloat)
        return new(data, resolution, resolution * size(data, 1), resolution * size(data, 2))
    end
end

"""
    # Parameters
    - src: ファイルパス
    - scan_size: スキャンサイズ (AA)
    - datatype: データ型
    - delimiter: 区切り文字
    - skipstart: 先頭何行を読み飛ばすか
    # Returns
    - Image
"""
function loadImage(
    src::String, scan_size::AbstractFloat;
    datatype::DataType=Float64, delimiter::Char='\t', skipstart::Integer=4
)::Image
    data = readdlm(src, delimiter, datatype; skipstart=skipstart)
    throw_if(!is_square(data), "Data must be square.")
    return Image(data; resolution=scan_size / size(data, 1))
end

"""
    # Parameters
    - from: 元画像
    - lowerleft: 左下の座標(ix,iy) (px)
    - size: 抽出するサイズ(x,y) (px)
    # Returns
    - Image
"""
function extract(from::Image, lowerleft::Tuple{Integer,Integer}, size::Tuple{Integer,Integer})::Image
    low = lowerleft[2]
    left = lowerleft[1]
    high = low + size[2] - 1
    right = left + size[1] - 1
    data = from.data[low:high, left:right]
    return Image(data; resolution = from.resolution)
end

"""
    # Summary
    - 与えられた行列でフィルタリングする
    # Parameters
    - from: 元画像
    - filter: フィルター 奇数x奇数の行列
    - trim: 端部を切り取るかどうか (default: true)
    # Returns
    - Image
"""
function filter(from::Image, filter::Matrix; trim::Bool=true)::Image
    if (size(filter, 1) % 2 == 0 || size(filter, 2) % 2 == 0)
        throw(ArgumentError("Filter size must be odd."))
    end
    halfheight, halfwidth = div.(size(filter), 2) # 1/2の切り捨て : halfX * 2 = X-1

    if trim
        size_ret = (size(from.data, 1) - 2 * halfheight, size(from.data, 2) - 2 * halfwidth)
    else
        size_ret = (size(from.data, 1), size(from.data, 2))
    end
    ret = Array{eltype(from.data)}(undef, size_ret...)

    for iy in halfheight+1:size(from.data, 1)-halfheight
        for ix in halfwidth+1:size(from.data, 2)-halfwidth
            window = from.data[iy-halfheight:iy+halfheight, ix-halfwidth:ix+halfwidth]
            ret[iy-iy_ini+1, ix-ix_ini+1] = sum(window .* filter)
        end
    end
    return Image(data, from.resolution)
end

"""
    # Summary
    - 横方向の微分を計算する
    # Parameters
    - from: 元画像
    # Returns
    - Image
"""
function horz_diff(from::Image)::Image
    filter = [0.5 -1 0.5]
    return filter(from, filter)
end

"""
    # Summary
    - 縦方向の微分を計算する
    # Parameters
    - from: 元画像
    # Returns
    - Image
"""
function vert_diff(from::Image)::Image
    filter = [0.5 -1 0.5]'
    return filter(from, filter)
end

########################################################
### Tip ################################################
########################################################

"""
    # Summary
    - チップデータを扱うための型
    # Fields
    - data: チップデータ
    - resolution: 解像度 (AA / px)
    - size: サイズ (AA)
    # Constructer Parameters
    - data: チップデータ
    - resolution: 解像度 (AA / px)
"""
struct Tip <: Surface
    data::Matrix
    resolution::AbstractFloat      # \AA / px
    size::AbstractFloat            # \AA : 正方形のみを許容

    function Tip(data::Matrix; resolution::AbstractFloat)::Tip
        throw_if(!is_square(data), "Data must be square.")
        return new(data, resolution, resolution * size(data, 1))
    end

end

end