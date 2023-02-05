"""
STMデータを扱うためのモジュール
"""
module SPMCore

using ..Unitful
using ..CUDA

export DEFAULT_UNIT
export Image, loadImage, extract, filter, horz_diff, vert_diff
export Tip

########################################################
### utilities ##########################################
########################################################

is_square(x::AbstractMatrix) = size(x, 1) == size(x, 2)

throw_if(x::Bool, msg::String) = x && throw(ArgumentError(msg))

const DEFAULT_UNIT = u"nm"

abstract type Surface end
(s::Surface)(x::Integer, y::Integer) = s.data[y, x]
(s::Surface)(xrange::UnitRange, yrange::UnitRange) = s.data[yrange, xrange]

########################################################
### Image ##############################################
########################################################

"""
    # Summary
    - 画像データを扱うための型
    # Fields
    - data: 画像データ (nm)
    - resolution: 解像度 (nm / px)
"""
struct Image <: Surface
    data::DenseMatrix{<:Real}
    resolution::AbstractFloat

    function Image(data::DenseMatrix{<:Real}, resolution::AbstractFloat; datatype=Float32)
        data = Matrix{datatype}(data)
        return new(data, resolution)
    end
end


"""
    # Parameters
    - from: 元画像
    - lowerleft: 左下の座標(ix,iy) (px)
    - size: 抽出するサイズ(x,y) (px)
    # Returns
    - Image
"""
function extract(from::Image, lowerleft::Tuple, extract_size::Tuple)::Image
    low = lowerleft[2]
    left = lowerleft[1]
    if (low < 1 || left < 1)
        throw(ArgumentError("Invalid lowerleft."))
    end

    high = low + extract_size[2] - 1
    right = left + extract_size[1] - 1
    if (high > size(from.data, 1) || right > size(from.data, 2))
        throw(ArgumentError("Invalid size."))
    end

    data = from.data[low:high, left:right]
    return Image(data, from.resolution)
end

"""
    # Summary
    - 1枚の画像から複数の同じサイズの画像を切り出す
    # Parameters
    - from: 元画像
    - lowerleft: 左下の座標(ix,iy) (px)を並べたベクトル
    - size: 抽出するサイズ(x,y) (px)
    # Returns
    - Vector{Image}
"""
function extract(
    from::Image, lowerleft::Vector, size::Tuple
)::Vector{Image}
    ret = Vector{Image}(undef, length(lowerleft))
    for i in 1:length(lowerleft)
        ret[i] = extract(from, lowerleft[i], size)
    end
    return ret
end

"""
    # Summary
    - 複数の画像から複数の同じサイズの画像を切り出す
    # Parameters
    - from: 元画像
    - lowerleft: 左下の座標(ix,iy) (px)のベクトルのベクトル
    - size: 抽出するサイズ(x,y) (px)
    - flatten: 1次元配列にするかどうか (default: false)
    # Returns
    - Vector{Image} or Vector{Vector{Image}}
"""
function extract(
    from::Vector{Image}, lowerleft::Vector, size::Tuple;
    flatten=false
)::Union{Vector{Image}, Vector{Vector{Image}}}
    if (length(lowerleft) != length(size))
        throw(ArgumentError("length(lowerleft) must be equal to length(size)."))
    end

    ret = Vector{Vector{Image}}(undef, length(from))
    for i in 1:length(from)
        ret[i] = extract(from[i], lowerleft[i], size)
    end

    if flatten
        return reduce(vcat, ret)
    else
        return ret
    end
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
function filter(from::Image, filter::DenseMatrix{<:Real}; trim::Bool=true)::Image
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
    - data: チップデータ (nm)
    - resolution: 解像度 (nm / px)
    # Constructor Parameters
    - data: チップデータ
    - resolution: 解像度 (nm / px)
"""
struct Tip <: Surface
    data::DenseMatrix{<:Real}
    resolution::AbstractFloat
end

end