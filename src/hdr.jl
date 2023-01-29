module HDR

using ..SPMCore

mutable struct HDRFile
    volume::Tuple{Integer, Integer}
    width::AbstractFloat
    height::AbstractFloat
    max_binval::Unsigned
    binval_bitlength::Unsigned
    max_realheight::AbstractFloat
end
HDRFile() = HDRFile((0, 0), 0.0, 0.0, 0, 0, 0.0)

function loadHDR(filename::String)::Image
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

    data = Matrix{datatype}(undef, file.volume...)
    dat_filename = replace(filename, ".hdr" => ".dat")
    open(dat_filename, "r") do io
        read!(io, data)
    end
    data = data * (file.max_realheight / file.max_binval) # nm

    resolution = file.width / size(data, 1) # nm / px

    return Image(data, resolution)
end

end