module HDR

using ..SPMCore
using Unitful

import Dates

# TODO : 上下がgwyddionと反転しているようなので修正したい

mutable struct HDRFile
    volume::Tuple{Integer, Integer}
    width::AbstractFloat
    height::AbstractFloat
    max_binval::Unsigned
    binval_bitlength::Unsigned
    max_realheight::AbstractFloat
end
HDRFile() = HDRFile((0, 0), 0.0, 0.0, 0, 0, 0.0)

function loadHDR(filename::String; outDataType=Float32)::Image
    if !occursin(".hdr", filename)
        throw(ArgumentError("File extension must be .hdr"))
    end
    file = HDRFile()
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
                state = :read_width
            elseif state == :read_width
                file.width = parse(Float64, split(line)[2])
                state = :wait_key
            elseif state == :wait_key && occursin("y_data", line)
                state = :wait_y_unit
            elseif state == :wait_y_unit
                state = :read_height
            elseif state == :read_height
                file.height = parse(Float64, split(line)[2])
                state = :wait_key
            elseif state == :wait_key && occursin("z_data", line)
                state = :wait_z_unit
            elseif state == :wait_z_unit
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
    outdata = Matrix{outDataType}(undef, size(raw_data')...)
    outdata .= raw_data' * (file.max_realheight / file.max_binval) # nm

    resolution = file.width / size(raw_data, 1) # nm / px

    return Image(outdata, resolution)
end

function saveHDR(
    surface::SPMCore.Surface, dist_path::String;
    unit::Unitful.FreeUnits=DEFAULT_UNIT,
    sample_name="", remark="", now=Dates.now(),
    overwrite=false, flat_threshold=1e-5
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
    if remark == ""
        remark = "Created by SPM.jl"
    end

    data = surface.data' .- minimum(surface.data)
    data[isnan.(data)] .= 0
    max_val = maximum(data)
    if max_val < flat_threshold
        max_converted = 1.0
        converted_data = zeros(UInt16, size(data))
    else
        max_converted = 2^16 - 1
        converted_data = round.(UInt16, data ./ max_val .* max_converted)
    end

    data_volume_y, data_volume_x = size(surface.data)
    y_end, x_end = size(surface.data) .* surface.resolution
    lines = [
        ":STM data",
        "\tformat vertion", "100",
        "\t:date; time",
            Dates.format(now, "yyyy/mm/dd"), Dates.format(now, "H:M:S"),
        "\t:sample name", sample_name,
        "\t:remark", remark,
        "\t:ascii flag; data type(2=byte; 3=+-byte; 4=word; 5=+-word; 8=float)", "0 4",
        "\t:data volume(x*y)", "$(data_volume_x)  $(data_volume_y)",
        "\t:x;y dimension(1=Length;2=Time;3=Current;4=Vt;5=Tmp;6=/L;7=/T;8=other)", "1 1",
        "\t:x_data -> unit; start; end; log flag", string(unit), "0 $(string(x_end)) 0",
        "\t:y_data -> unit; start; end; log flag", string(unit), "0 $(string(y_end)) 0",
        "\t:z_data -> unit; max; min; (max; mini)(before conversion); log flag",
            string(unit), "$(max_converted) 0 $(max_val) 0.0 0",
        "\t:reserved", "0 0 0 0",
        "\t:stm voltage unit; stm current unit; A/D name", "V", "nA", "unknown",
        "a"
    ]
    open(dist_path, "w") do f
        for line in lines
            println(f, line*"\r")
        end
    end

    dat_path = replace(dist_path, ".hdr" => ".dat")
    open(dat_path, "w") do f
        write(f, converted_data)
    end

    return true
end

end