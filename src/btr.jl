module BTR

using ..SPMCore
import ..HDR
import Statistics.mean
import MDToolbox as MDT
import Flux
import Dates

# TODO : chainrulescore.rruleが定義されているか確認

########################################################
### BTR ################################################
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
    px_size::Integer, resolution::AbstractFloat;
    datatype::DataType=Float64
)::Tip
    return Tip(zeros(datatype, px_size, px_size), resolution)
end

function genFlatTip(px_size::Integer, image::Image)::Tip
    return genFlatTip(px_size, image.resolution; datatype=eltype(image.data))
end

abstract type BTRResult end

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
        threshold = replace(string(result.threshold), "." => "_")
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
    loss_history::Vector{Real}
end

function saveResult(
    result::DifferentiableBTRResult, filename::String; timestamp::Bool=true, savehistory::Bool=true
)::Bool
    if timestamp
        filename = filename * "_" * Dates.format(Dates.now(), "yyyymmdd_HMS")
    end

    lambda = replace(string(result.lambda), "." => "_")
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
        lambda = replace(string(result.lambda), "." => "_")
        saveResult(result, joinpath(dirname, "result_lambda_$(lambda)"), timestamp=false, savehistory=false)
    end

    max_epoch = results[1].max_epoch
    open(joinpath(dirname, "loss_history.csv"), "w") do io
        header = "epoch, " * join(["lambda=$(replace(string(result.lambda), "." => "_"))" for result in results], ", ")
        println(io, header)
        for i in 1:max_epoch
            loss = join([string(result.loss_history[i]) for result in results], ", ")
            println(io, "$(i), $(loss)")
        end
    end

    return true
end

Flux.@functor Tip (data,)
(m::Tip)(image_data) = MDT.idilation(MDT.ierosion(image_data, m.data), m.data)

function solveDifferentiableBTR(
        images::Vector{Image}, tip_size::Integer, max_epoch::Integer, lambda::Real;
        debug_interval=20
)::DifferentiableBTRResult
    images_copy = deepcopy(images)

    tip = genFlatTip(tip_size, images[1].resolution)
    min_loss_tip = deepcopy(tip)

    loss(vim_cp::Vector{Image}, vim::Vector{Image}) = mean(Flux.Losses.mse.(tip.([im.data for im in vim_cp]), [im.data for im in vim]))

    # 訓練するモデルを与える
    ps = Flux.params(tip)

    # 訓練データを与える
    train_loader = Flux.Data.DataLoader(
        (data=images_copy, label=images), batchsize=1, shuffle=false
    )

    # 最適化手法を指定する
    opt = Flux.ADAMW(1.0, (0.9, 0.999), lambda)
    @info "$(Threads.threadid())th thread : optimizer setup completed"

    loss_train = Vector{Float64}(undef, max_epoch)
    min_loss = Inf
    t_ini = Dates.now()
    for epoch in 1:max_epoch
        for (x, y) in train_loader
            gs = Flux.gradient(() -> loss(x, y), ps)
            Flux.Optimise.update!(opt, ps, gs)
            tip.data .= min.(tip.data, 0.0)
            tip.data .= MDT.translate_tip_mean(tip.data)
        end
        epoch_loss = loss(images_copy, images)
        if epoch_loss < min_loss
            min_loss = epoch_loss
            min_loss_tip = deepcopy(tip)
        end

        loss_train[epoch] = loss(images_copy, images)

        if epoch % debug_interval == 0
            @info "$(Threads.threadid())th thread : $(epoch)th epoch completed in $((Dates.now() - t_ini).value/1000) sec"
        end
    end

    DifferentiableBTRResult(lambda, max_epoch, tip, min_loss_tip, loss_train)
end

function solveDifferentiableBTR(
    images::Vector{Image}, tip_size::Integer, max_epoch::Integer, lambdas::Vector;
    debug_interval=20
)::Vector{DifferentiableBTRResult}
    ret = Vector{DifferentiableBTRResult}(undef, length(lambdas))
    Threads.@threads for it in eachindex(lambdas)
        if typeof(lambdas[it]) <: Real
            @info "$(Threads.threadid())th thread : start solving for lambda = $(lambdas[it])"
            ret[it] = solveDifferentiableBTR(images, tip_size, max_epoch, lambdas[it]; debug_interval=debug_interval)
        else
            @warn "lambdas[$(it)] is not a Real. skipping..."
            continue
        end
        @info "$(Threads.threadid())th thread : $(it)th lambda completed"
    end
    return ret
end

end