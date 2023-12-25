module HDR

using ..SPMCore

import Dates

export loadHDR, saveHDR

mutable struct HDRFile
    volume::Tuple{Integer, Integer}
    width::AbstractFloat
    height::AbstractFloat
    max_binval::Unsigned
    binval_bitlength::Unsigned
    max_realheight::AbstractFloat
    horz_unit::String
    vert_unit::String
end
HDRFile() = HDRFile((0, 0), 0.0, 0.0, 0, 0, 0.0, "nm", "nm")


"""
    Parameters
    - filename: ファイル名
    Returns
    - surface: 表面データ
"""
function loadHDR(filename::String)::SPMCore.Surface
    if !occursin(".hdr", filename)
        throw(ArgumentError("File extension must be .hdr"))
    end
    file = HDRFile()
    x_unit, y_unit, z_unit = "", "", ""
    open(filename, "r") do f
        state = :wait_key
        for line in eachline(f)
            if state == :wait_key && occursin("volume(x*y)", line)
                state = :read_volume
            elseif state == :read_volume
                vals = split(line)
                file.volume = (parse(Int, vals[1]), parse(Int, vals[2]))
                state = :wait_key
            elseif state == :wait_key && occursin("x_data", line)
                state = :wait_x_unit
            elseif state == :wait_x_unit
                x_unit = split(line)[1]
                state = :read_width
            elseif state == :read_width
                file.width = parse(Float64, split(line)[2])
                state = :wait_key
            elseif state == :wait_key && occursin("y_data", line)
                state = :wait_y_unit
            elseif state == :wait_y_unit
                y_unit = split(line)[1]
                state = :read_height
            elseif state == :read_height
                file.height = parse(Float64, split(line)[2])
                state = :wait_key
            elseif state == :wait_key && occursin("z_data", line)
                state = :wait_z_unit
            elseif state == :wait_z_unit
                z_unit = split(line)[1]
                state = :read_z_data
            elseif state == :read_z_data
                splitted = split(line)
                file.max_binval = parse(Int, splitted[1])
                file.binval_bitlength = convert(Unsigned, log2(file.max_binval+1))
                file.max_realheight = parse(Float64, splitted[3])
                state = :exit
            elseif state == :exit
                break
            end
        end
    end

    if x_unit != y_unit
        throw(ArgumentError("Inconsistent unit: $x_unit, $y_unit"))
    elseif !in(x_unit, SPMCore.SUPPORTED_UNITS)
        throw(ArgumentError("Unsupported unit: $x_unit"))
    end
    file.horz_unit = x_unit

    if !in(z_unit, SPMCore.SUPPORTED_UNITS)
        throw(ArgumentError("Unsupported unit: $z_unit"))
    end
    file.vert_unit = z_unit

    if file.binval_bitlength == 8
        datatype = UInt8
    elseif file.binval_bitlength == 16
        datatype = UInt16
    else
        throw(ArgumentError("Unsupported bit length: $(file.binval_bitlength)"))
    end

    raw_data = Matrix{datatype}(undef, file.volume...)
    dat_filename = replace(filename, ".hdr" => ".dat")
    open(dat_filename, "r") do io
        read!(io, raw_data)
    end
    outdata = Matrix{Float32}(undef, size(raw_data')...)
    outdata .= raw_data' * (file.max_realheight / file.max_binval) # nm

    resolution = file.width / size(raw_data, 1) # nm / px

    return Surface(outdata, file.vert_unit, Float32(resolution), file.horz_unit)
end

"""
    Parameters
    - surface: 表面データ
    - dist_path: 保存先のパス
    - sample_name: サンプル名
    - remark: 備考
    - now: 作成日時
    - overwrite: 上書きするかどうか
    - flat_threshold: 平坦化する閾値
    Returns
    - true: 成功
"""
function saveHDR(
    surface::SPMCore.Surface, dist_path::String;
    sample_name="", remark="", now=Dates.now(),
    overwrite=false
)::Bool
    if !occursin(".hdr", dist_path)
        throw(ArgumentError("File extension must be .hdr"))
    elseif isfile(dist_path)
        if overwrite
            @info "File already exists: $(dist_path). Overwrite."
        else
            throw(ArgumentError("File already exists: $(dist_path)"))
        end
    end

    # ヘッダファイルの書き込み
    if sample_name == ""
        sample_name = replace(basename(dist_path), ".hdr" => "")
    end
    remark = remark * " - Created by SPM.jl"

    data = surface.matrix' .- minimum(surface.matrix)
    data[isnan.(data)] .= 0
    # 単位はnmに統一
    out_unit = "nm"
    if surface.vert_unit == "pm"
        data = data * 1e-3
    elseif surface.vert_unit == "AA"
        data = data * 1e-1
    elseif surface.vert_unit == "um"
        data = data * 1e3
    end

    max_val = maximum(data)
    max_converted = 2^16 - 1
    converted_data = round.(UInt16, data ./ max_val .* max_converted)

    data_volume_y, data_volume_x = size(surface.matrix)
    horz_resolution = surface.horz_resolution
    if surface.horz_unit == "pm"
        horz_resolution = horz_resolution * 1e-3
    elseif surface.horz_unit == "AA"
        horz_resolution = horz_resolution * 1e-1
    elseif surface.horz_unit == "um"
        horz_resolution = horz_resolution * 1e3
    end
    y_end, x_end = size(surface.matrix) .* horz_resolution
    lines = [
        ":STM data",
        "\tformat vertion",
            "100",
        "\t:date; time",
            Dates.format(now, "yyyy/mm/dd"),
            Dates.format(now, "H:M:S"),
        "\t:sample name",
            sample_name,
        "\t:remark",
            remark,
        "\t:ascii flag; data type(2=byte; 3=+-byte; 4=word; 5=+-word; 8=float)",
            "0 4",
        "\t:data volume(x*y)",
            "$(data_volume_x)  $(data_volume_y)",
        "\t:x;y dimension(1=Length;2=Time;3=Current;4=Vt;5=Tmp;6=/L;7=/T;8=other)",
            "1 1",
        "\t:x_data -> unit; start; end; log flag",
            string(out_unit),
            "0 $(string(x_end)) 0",
        "\t:y_data -> unit; start; end; log flag",
            string(out_unit),
            "0 $(string(y_end)) 0",
        "\t:z_data -> unit; max; min; (max; mini)(before conversion); log flag",
            string(out_unit),
            "$(max_converted) 0 $(max_val) 0.0 0",
        "\t:reserved",
            "0 0 0 0",
        "\t:stm voltage unit; stm current unit; A/D name",
            "V",
            "nA",
            "unknown",
        "a"
    ]

    open(dist_path, "w") do f
        for line in lines
            println(f, line*"\r")
        end
    end

    dat_path = replace(dist_path, ".hdr" => ".dat")
    open(dat_path, "w") do f
        write(f, converted_data) # バイナリ書き込み
    end

    return true
end

end # module HDR