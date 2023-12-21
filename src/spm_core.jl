"""
STMデータを扱うためのモジュール
"""
module SPMCore

export SUPPORTED_UNITS, DEFAULT_UNIT
export Surface, extract, convolute

########################################################
### utilities ##########################################
########################################################

throw_if(x::Bool, msg::String) = x && throw(ArgumentError(msg))

SUPPORTED_UNITS = ("pm", "AA", "nm", "um")

const DEFAULT_UNIT = SUPPORTED_UNITS[3] # nm

########################################################
### Surface ############################################
########################################################

"""
    Summary
    - 2次元の表面データ
    Parameters
    - matrix: 表面データ
    - vert_unit: 垂直方向の単位
    - horz_resolution: 水平方向の解像度
    - horz_unit: 水平方向の単位
"""
struct Surface
    matrix::Matrix{Float32}    # vert_unit で表した高さの配列
    vert_unit::String    # 垂直方向の単位
    horz_resolution::Float32        # horz_unit / px で表した水平方向の解像度
    horz_unit::String    # 水平方向の単位

    function Surface(
        matrix::Matrix{Float32}, vert_unit::String, horz_resolution::Float32, horz_unit::String
    )
        throw_if(!(vert_unit in SUPPORTED_UNITS), "Unsupported unit: $vert_unit")
        throw_if(!(horz_unit in SUPPORTED_UNITS), "Unsupported unit: $horz_unit")
        return new(matrix, vert_unit, horz_resolution, horz_unit)
    end

    function Surface(matrix::Matrix{Float32}, resolution::Float32)
        return Surface(matrix, DEFAULT_UNIT, resolution, DEFAULT_UNIT)
    end
end

(s::Surface)(x::Integer, y::Integer) = s.matrix[y, x]   # surface(x,y) でs.matrix[y,x]を返す
(s::Surface)(x_range::UnitRange, y_range::UnitRange) = s.matrix[y_range, x_range]


"""
    Parameters
    - from: 元画像
    - lowerleft: 左下の座標(ix,iy) (px)
    - size: 抽出するサイズ(x,y) (px)
    Returns
    - Image
"""
function extract(from::Surface, lowerleft::Tuple, extract_size::Tuple)::Surface
    low = lowerleft[2]
    left = lowerleft[1]
    if (low < 1 || left < 1)
        throw(ArgumentError("Invalid lowerleft."))
    end

    high = low + extract_size[2] - 1
    right = left + extract_size[1] - 1
    if (high > size(from.matrix, 1) || right > size(from.matrix, 2))
        throw(ArgumentError("Invalid size."))
    end

    matrix = from(low:high, left:right)
    return Surface(matrix, from.vert_unit, from.horz_resolution, from.horz_unit)
end


"""
    Summary
    - 与えられた行列で畳み込みを行う
    Parameters
    - from: 元画像
    - filter: フィルター 奇数x奇数の行列
    Returns
    - Image
"""
function convolute(from::Surface, filter::Matrix)::Surface
    if (size(filter, 1) % 2 == 0 || size(filter, 2) % 2 == 0)
        throw(ArgumentError("Filter size must be odd."))
    end
    halfheight, halfwidth = div.(size(filter), 2) # 1/2の切り捨て : halfX * 2 = X-1

    # 端部ではフィルターのサイズが足りないので鏡像で埋める
    mat = Matrix{eltype(from.matrix)}(undef, size(from.matrix, 1)+2*halfheight, size(from.matrix, 2)+2*halfwidth)
    mat[halfheight+1:end-halfheight, halfwidth+1:end-halfwidth] .= from.matrix
    mat[:, 1:halfwidth] .= mat[:, 2*halfwidth+1:-1:halfwidth+2]
    mat[:, end-halfwidth+1:end] .= mat[:, end-halfwidth-2:-1:end-2*halfwidth-1]
    mat[1:halfheight, :] .= mat[2*halfheight+1:-1:halfheight+2, :]
    mat[end-halfheight+1:end, :] .= mat[end-halfheight-2:-1:end-2*halfheight-1, :]

    ret = Matrix{eltype(from.matrix)}(undef, size(from.matrix)...)

    for iy in halfheight+1:halfheight+size(from.matrix, 1)
        for ix in halfwidth+1:halfwidth+size(from.matrix, 2)
            ret[iy-halfheight, ix-halfwidth] = sum(mat[iy-halfheight:iy+halfheight, ix-halfwidth:ix+halfwidth] .* filter)
        end
    end
    return Surface(ret, from.vert_unit, from.horz_resolution, from.horz_unit)
end

end # module SPMCore