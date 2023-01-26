module BTR

using ..SPMCore
import Statistics.mean
import MDToolbox as MDT
import Flux

########################################################
### BTR ################################################
########################################################

match(image::Image, tip::Tip) = image.resolution == tip.resolution

function ierosion(image::Image, tip::Tip)::Image
    if !match(image, tip)
        throw(ArgumentError("image and tip must have the same resolution"))
    end
    return Image(MDT.ierosion(image.data, tip.data); resolution=image.resolution)
end

function idilation(image::Image, tip::Tip)::Image
    if !match(image, tip)
        throw(ArgumentError("image and tip must have the same resolution"))
    end
    return Image(MDT.idilation(image.data, tip.data); resolution=image.resolution)
end

function opening(image::Image, tip::Tip)::Image
    return idilation(ierosion(image, tip), tip)
end

function genFlatTip(
        px_size::Integer, resolution::AbstractFloat;
        datatype::DataType=Float64
)::Tip
    return Tip(zeros(datatype, px_size, px_size); resolution=resolution)
end

function genFlatTip(px_size::Integer, image::Image)::Tip
    return genFlatTip(px_size, image.resolution; datatype=eltype(image.data))
end

abstract type BTRResult end

struct VillarrubiaBTRReslt <: BTRResult
    tip::Tip
    loss::AbstractFloat
    threshold::AbstractFloat
end
OriginalBTRResult = VillarrubiaBTRReslt

function solveVillarrubiaBTR(
        images::Vector{Image}, tip_size::Integer, threshold::AbstractFloat
)::OriginalBTRResult
    P = zeros(eltype(images[1].data), tip_size, tip_size)
    imageData = [image.data for image in images]
    MDT.itip_estimate!(P, imageData; thresh=threshold)
    tip = Tip(P; resolution=images[1].resolution)

    loss = 0.0
    for image in images
        loss += mean((opening(image, tip).data .- image.data) .^ 2)
    end
    return VillarrubiaBTRReslt(tip, loss, threshold)
end
solveOriginalBTR = solveVillarrubiaBTR

function solveVillarrubiaBTRoverThresholds(
        images::Vector{Image}, tip_size::Integer, thresholds::Vector{AbstractFloat}
)::Vector{OriginalBTRResult}
    ret = Vector{OriginalBTRResult}(undef, length(thresholds))
    Threads.@threads for it in eachindex(thresholds)
        ret[it] = solveVillarrubiaBTR(images, tip_size, thresholds[it])
    end
    return ret
end
solveOriginalBTRoverThresholds = solveVillarrubiaBTRoverThresholds

struct DifferentiableBTRResult <: BTRResult
    lambda::AbstractFloat
    max_epoch::Integer
    min_loss::AbstractFloat
    final_tip::Tip
    loss_minimizing_tip::Tip
end

(m::Tip)(image::Image)::Image = opening(image, m)

function solveDifferentiableBTR(
        images::Vector{Image}, tip_size::Integer,
        max_epoch::Integer, lambda::AbstractFloat
)
    images_copy = deepcopy(images)

    tip = genFlatTip(tip_size, images[1].resolution)
    min_loss_tip = deepcopy(tip)

    Flux.@functor Tip (data,)
    function loss(images_copy::Vector{Image}, images::Vector{Image})
        data_opened = [tip(image).data for image in images_copy]
        data = [image.data for image in images]
        squared_error = [ mean((data_opened[it] .- data[it]).^2) for it in eachindex(data) ]
        return mean(squared_error)
    end

    # 訓練するモデルを与える
    ps = Flux.params(tip)

    # 訓練データを与える
    train_loader = Flux.Data.DataLoader(
        (data=images_copy, label=images), batchsize=1, shuffle=false
    )

    # 最適化手法を指定する
    opt = Flux.ADAMW(1.0, (0.9, 0.999), lambda)

    loss_train = Vector{Float64}(undef, max_epoch)
    min_loss = Inf
    for epoch in 1:max_epoch
        for (x, y) in train_loader
            @info typeof(x), typeof(y)
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
    end

    DifferentiableBTRResult(lambda, max_epoch, min_loss, tip, min_loss_tip)
end

end