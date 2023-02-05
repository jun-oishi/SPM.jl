module BTR

using ..SPMCore
using ..CUDA

import ..HDR
import Statistics.mean
import MDToolbox as MDT
import Flux
import Dates
import Printf.@sprintf

# TODO : chainrulescore.rruleが定義されているか確認

########################################################
### utilities ##########################################
########################################################

match(image::Image, tip::Tip) = image.resolution == tip.resolution

function ierosion(image::Image, tip::Tip)::Image
    if !match(image, tip)
        throw(ArgumentError("image and tip must have the same resolution"))
    end
    return Image(MDT.ierosion(image.data, tip.data), image.resolution)
end

function idilation(image::Image, tip::Tip)::Image
    if !match(image, tip)
        throw(ArgumentError("image and tip must have the same resolution"))
    end
    return Image(MDT.idilation(image.data, tip.data), image.resolution)
end

function opening(image::Image, tip::Tip)::Image
    return idilation(ierosion(image, tip), tip)
end

function genFlatTip(
    px_size::Integer, resolution::Real;
    datatype::DataType=Float64, use_cuda::Bool=false
)::Tip
    if use_cuda
        data = CUDA.zeros(datatype, px_size, px_size)
    else
        data = zeros(datatype, px_size, px_size)
    end
    return Tip(data, resolution)
end

function genFlatTip(px_size::Integer, image::Image)::Tip
    return genFlatTip(px_size, image.resolution; datatype=eltype(image.data))
end

"""
    # Summary
    - 行列をratio倍に拡大する
    # Parameters
    - mat: 拡大する行列
    - ratio: 拡大倍率
    # Returns
    - matの行数と列数をratio倍に拡大した行列
"""
function upscaleMatrix(mat::AbstractMatrix, ratio::Integer)::AbstractMatrix
    return [ mat[div(i,ratio)+1, div(j,ratio)+1] for i in 0:ratio*size(mat, 1)-1, j in 0:ratio*size(mat, 2)-1 ]
end

abstract type BTRResult end

########################################################
### Villarrubia ########################################
########################################################

struct VillarrubiaBTRResult <: BTRResult
    tip::Tip
    loss::AbstractFloat
    threshold::AbstractFloat
end
OriginalBTRResult = VillarrubiaBTRResult

function saveResult(
    result::VillarrubiaBTRResult, filename::String; timestamp::Bool=true
)::Bool
    sample_name = filename
    if timestamp
        filename = filename * "_" * Dates.format(Dates.now(), "yyyymmdd_HMS")
        sample_name = filename
    end
    filename = filename * ".hdr"

    remark = "Villarrubia BRT result: threshold=$(result.threshold) => loss=$(result.loss)"
    HDR.saveHDR(result.tip, filename; sample_name = sample_name, remark=remark)
    return true
end

function saveResults(
    source_images::Vector{Image}, results::Vector{VillarrubiaBTRResult}, dirname::String; timestamp::Bool=true
)::Bool
    if timestamp
        dirname = dirname * "_" * Dates.format(Dates.now(), "yyyymmdd_HMS")
    end
    mkdir(dirname)

    for (i, image) in enumerate(source_images)
        HDR.saveHDR(image, joinpath(dirname, "source_$(i).hdr"))
    end

    for result in results
        threshold = replace(@sprintf("%.3e", result.threshold), "." => "_")
        saveResult(result, joinpath(dirname, "result_threshold_$(threshold)"), timestamp=false)
    end
    return true
end

function solveVillarrubiaBTR(
    images::Vector{Image}, tip_size::Integer, threshold::Real
)::OriginalBTRResult
    P = zeros(eltype(images[1].data), tip_size, tip_size)
    imageData = [image.data for image in images]
    MDT.itip_estimate!(P, imageData; thresh=threshold)
    tip = Tip(P, images[1].resolution)

    loss = 0.0
    for image in images
        loss += mean((opening(image, tip).data .- image.data) .^ 2)
    end
    return VillarrubiaBTRResult(tip, loss, threshold)
end

function solveVillarrubiaBTR(
    images::Vector{Image}, tip_size::Integer, thresholds::Vector
)::Vector{VillarrubiaBTRResult}
    ret = Vector{VillarrubiaBTRResult}(undef, length(thresholds))
    Threads.@threads for it in eachindex(thresholds)
        if typeof(thresholds[it]) <: Real
            @info "$(Threads.threadid())th thread : start solving for thresh = $(thresholds[it])"
            ret[it] = solveVillarrubiaBTR(images, tip_size, thresholds[it])
        else
            @warn "thresholds[$(it)] is not a Real number. skipping..."
            continue
        end
    end
    return ret
end
solveOriginalBTR = solveVillarrubiaBTR

struct DifferentiableBTRResult <: BTRResult
    lambda::Real
    max_epoch::Integer
    tip::Tip        # 最終形状
    loss_minimizing_tip::Tip
    loss_history::Vector{<:Real}
end

function saveResult(
    result::DifferentiableBTRResult, filename::String; timestamp::Bool=true, savehistory::Bool=true
)::Bool
    if timestamp
        filename = filename * "_" * Dates.format(Dates.now(), "yyyymmdd_HMS")
    end

    lambda = @sprintf("%.8f", result.lambda)
    max_epoch = result.max_epoch
    loss_final = replace(string(result.loss_history[end]), "." => "_")
    remark = "result of differentiable BTR: lambda=$(lambda), epoch=$(max_epoch) => loss=$(loss_final)"
    HDR.saveHDR(result.tip, filename*"_final.hdr"; remark=remark)

    min_loss_epoch = argmin(result.loss_history)
    min_loss = replace(string(result.loss_history[min_loss_epoch]), "." => "_")
    remark = "result of differentiable BTR: lambda=$(lambda), epoch=$(min_loss_epoch) => loss=$(min_loss)"
    HDR.saveHDR(result.loss_minimizing_tip, filename*"_min_loss.hdr"; remark=remark)

    if savehistory
        unit = string(DEFAULT_UNIT)
        open(filename*"_loss_history.csv", "w") do io
            println(io, "epoch, loss[$(unit)^2]")
            for (i, loss) in enumerate(result.loss_history)
                println(io, "$(i), $(loss)")
            end
        end
    end

    return true
end

function saveResults(
    source_images::Vector{Image}, results::Vector{DifferentiableBTRResult}, dirname::String; timestamp::Bool=true
)::Bool
    if timestamp
        dirname = dirname * "_" * Dates.format(Dates.now(), "yyyymmdd_HMS")
    end
    mkdir(dirname)

    for (i, image) in enumerate(source_images)
        HDR.saveHDR(image, joinpath(dirname, "source_$(i).hdr"))
    end

    for result in results
        lambda = replace(@sprintf("%.3e", result.lambda), "." => "_")
        saveResult(result, joinpath(dirname, "result_lambda_$(lambda)"), timestamp=false, savehistory=false)
    end

    max_epoch = results[1].max_epoch
    open(joinpath(dirname, "loss_history.csv"), "w") do io
        header = "epoch, " * join(["lambda=$(@sprintf("%.8e", result.lambda))" for result in results], ", ")
        println(io, header)
        for i in 1:max_epoch
            loss = join([string(result.loss_history[i]) for result in results], ", ")
            println(io, "$(i), $(loss)")
        end
    end

    return true
end

Flux.@functor Tip (data,)

function solveDifferentiableBTR(
    images::Vector{Image}, tip_size::Integer, max_epoch::Integer, lambda::Real;
    debug_interval=20, need_loss_minimizing=false, use_cuda=false, downscale_ratio::Integer=1
)::DifferentiableBTRResult
    # TODO downscaleの実装
    image_resolution = images[1].resolution
    for image in images
        if image.resolution != image_resolution
            @error "resolution of images are not same. aborting..."
            return
        end
    end

    tip_resolution = image_resolution * downscale_ratio

    if use_cuda
        image_data = Vector{CuMatrix{Float32}}([Float32.(image.data) for image in images])
        image_data_copy = deepcopy(image_data)

        tip = genFlatTip(tip_size, tip_resolution; datatype=Float32, use_cuda=true)
    else
        d_type = typeof(images[1].data[1, 1])
        image_data = Vector{Matrix{d_type}}([image.data for image in images])
        image_data_copy = deepcopy(image_data)

        tip = genFlatTip(tip_size, tip_resolution; datatype=d_type, use_cuda=false)
    end

    function loss(v_mat_cp, v_mat)
        tip_data = upscaleMatrix(tip.data, downscale_ratio)
        opened = [MDT.idilation(MDT.ierosion(mat_cp, tip_data), tip_data) for mat_cp in v_mat_cp]
        return mean(Flux.Losses.mse.(opened, v_mat))
    end

    # 訓練するモデルを与える
    ps = Flux.params(tip)

    # 訓練データを与える
    train_loader = Flux.Data.DataLoader(
        (data=image_data_copy, label=image_data), batchsize=1, shuffle=false
    )

    # 最適化手法を指定する
    opt = Flux.ADAMW(1.0, (0.9, 0.999), lambda)
    @info "$(Threads.threadid())th thread : optimizer setup completed"

    debug_interval = debug_interval == 0 ? max_epoch : debug_interval
    n_interval = div(max_epoch, debug_interval)
    max_epoch = n_interval * debug_interval
    loss_train = Vector{Float64}(undef, max_epoch)
    min_loss_tip = deepcopy(tip)
    min_loss = Inf
    epoch = 0

    t_ini = Dates.now()
    if need_loss_minimizing
        for interval in 1:n_interval
            for i in 1:debug_interval
                epoch += 1
                for (x, y) in train_loader
                    gs = Flux.gradient(() -> loss(x, y), ps)
                    Flux.Optimise.update!(opt, ps, gs)
                    tip.data .= min.(tip.data, 0.0)
                    tip.data .= MDT.translate_tip_mean(tip.data)
                end
                loss_train[epoch] = loss(image_data_copy, image_data)
                if loss_train[epoch] < min_loss
                    min_loss = loss_train[epoch]
                    min_loss_tip.data .= tip.data
                end
            end
            @info "$(Threads.threadid())th thread : $(epoch)th epoch completed in $((Dates.now() - t_ini).value/1000) sec"
        end
    else
        for interval in 1:n_interval
            for i in 1:debug_interval
                epoch += 1
                for (x, y) in train_loader
                    gs = Flux.gradient(() -> loss(x, y), ps)
                    Flux.Optimise.update!(opt, ps, gs)
                    tip.data .= min.(tip.data, 0.0)
                    tip.data .= MDT.translate_tip_mean(tip.data)
                end
                loss_train[epoch] = loss(image_data_copy, image_data)
            end
            @info "$(Threads.threadid())th thread : $(epoch)th epoch completed in $((Dates.now() - t_ini).value/1000) sec"
        end
    end

    DifferentiableBTRResult(lambda, max_epoch, tip, min_loss_tip, loss_train)
end

function solveDifferentiableBTR(
    images::Vector{Image}, tip_size::Integer, max_epoch::Integer, lambdas::Vector{<:Real};
    save_on_each_lambda=false, kw...
)::Vector{DifferentiableBTRResult}
    ret = Vector{DifferentiableBTRResult}(undef, length(lambdas))
    Threads.@threads for it in eachindex(lambdas)
        @info "$(Threads.threadid())th thread : start solving for lambda = $(lambdas[it])"
        ret[it] = solveDifferentiableBTR(
            images, tip_size, max_epoch, lambdas[it]; kw...
        )
        if save_on_each_lambda
            dirname = "result_$(replace(@sprintf("%.3e", lambdas[it]), "." => "_"))"
            if !isdir(dirname)
                mkpath(dirname)
            end
            saveResult(dirname, ret[it])
        end
        @info "$(Threads.threadid())th thread : $(it)th lambda completed"
    end
    return ret
end

end