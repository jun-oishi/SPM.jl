module BTR

using ..SPMCore
# using ..CUDA

import ..HDR
import Statistics
import MDToolbox as MDT
import Flux
import Dates
import Printf.@sprintf

# TODO : chainrulescore.rruleが定義されているか確認

########################################################
### utilities ##########################################
########################################################

ierosion(image::Surface, tip::Surface)::Surface = Surface(MDT.ierosion(image.matrix, tip.matrix), mold=image)
idilation(image::Surface, tip::Surface)::Surface = Surface(MDT.idilation(image.matrix, tip.matrix), mold=image)
opening(image::Surface, tip::Surface)::Surface = idilation(ierosion(image, tip), tip)
refineImage(image::Surface, tip::Surface)::Surface = ierosion(image, tip)

abstract type BTRResult end


########################################################
### Villarrubia ########################################
########################################################

"""
    Summary
    - VillarrubiaのBTRの結果
"""
struct VillarrubiaBTRResult <: BTRResult
    tip::Surface
    loss::AbstractFloat
    threshold::AbstractFloat
end

"""
    Summary
    - VillarrubiaのBTRの結果を保存する
    Parameters
    - result: 結果
    - filename: 保存するファイル名
    Returns
    - Bool: 正常終了したらtrue
"""
function saveResult(
    result::VillarrubiaBTRResult, filename::String
)::Bool
    sample_name = filename
    filename = filename * "_" * Dates.format(Dates.now(), "yyyymmdd_HMS")
    sample_name = filename
    filename = filename * ".hdr"

    remark = "Villarrubia BRT result: threshold=$(result.threshold) => loss=$(result.loss)"
    HDR.saveHDR(result.tip, filename; sample_name = sample_name, remark=remark)
    return true
end

"""
    Summary
    - VillarrubiaのBTRを解いて探針形状を求める
    Parameters
    - images: 元画像の配列
    - tip_size: 探針のサイズ[px]
    - threshold: 閾値
    Returns
    - VillarrubiaBTRResult
"""
function solveVillarrubiaBTR(
    images::Vector{Surface}, tip_size::Integer, threshold::Real
)::VillarrubiaBTRResult
    P = zeros(Float32, tip_size, tip_size)
    imageData = [image.matrix for image in images]
    MDT.itip_estimate!(P, imageData; thresh=threshold)
    tip = Surface(P, mold=images[1])

    loss = 0.0
    for image in images
        loss += Statistics.mean((opening(image, tip).matrix .- image.matrix) .^ 2)
    end
    return VillarrubiaBTRResult(tip, loss, threshold)
end


########################################################
### Differentiable BTR #################################
########################################################

"""
    Summary
    - Differentiable BTRの結果
"""
struct DifferentiableBTRResult <: BTRResult
    lambda::Real
    max_epoch::Integer
    tip::Surface        # 最終形状
    loss_minimizing_tip::Surface
    loss_history::Vector{<:Real}
end

"""
    Summary
    - Differentiable BTRの結果を保存する
        最終形状、最小loss時の形状(定義されていれば)、lossの履歴を保存する
    Parameters
    - result: 結果
    - filename: 保存するファイル名
    - timestamp: タイムスタンプをつけるかどうか
    - savehistory: trueならlossの履歴をcsvで保存する
    Returns
    - Bool: 正常終了したらtrue
"""
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

    if !isempty(result.loss_minimizing_tip.matrix)
        min_loss_epoch = argmin(result.loss_history)
        min_loss = replace(string(result.loss_history[min_loss_epoch]), "." => "_")
        remark = "result of differentiable BTR: lambda=$(lambda), epoch=$(min_loss_epoch) => loss=$(min_loss)"
        HDR.saveHDR(result.loss_minimizing_tip, filename*"_min_loss.hdr"; remark=remark)
    end

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

Flux.@functor Surface (matrix,)

function solveDifferentiableBTR(
    images::Vector{Surface}, tip_size::Integer, max_epoch::Integer, lambda::Real;
    debug_interval=20, need_loss_minimizing=false
)::DifferentiableBTRResult

    image_data = Vector{Matrix{Float32}}([image.matrix for image in images])
    image_data_copy = deepcopy(image_data)

    tip = Surface(zeros(Float32, tip_size, tip_size), mold=images[1])

    function loss(v_mat_cp, v_mat)
        tip_data = tip.matrix
        opened = [MDT.idilation(MDT.ierosion(mat_cp, tip_data), tip_data) for mat_cp in v_mat_cp]
        return Statistics.mean(Flux.Losses.mse.(opened, v_mat))
    end

    # 訓練するモデルを与える
    ps = Flux.params(tip)

    # 訓練データを与える
    train_loader = Flux.Data.DataLoader(
        (data=image_data_copy, label=image_data), batchsize=1, shuffle=false
    )

    # 最適化手法を指定する
    opt = Flux.ADAMW(1.0, (0.9, 0.999), lambda)
    # @info "$(Threads.threadid())th thread : optimizer setup completed"

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
                    tip.matrix .= min.(tip.matrix, 0.0)
                    tip.matrix .= MDT.translate_tip_mean(tip.matrix)
                end
                loss_train[epoch] = loss(image_data_copy, image_data)
                if loss_train[epoch] < min_loss
                    min_loss = loss_train[epoch]
                    min_loss_tip.matrix .= tip.matrix
                end
            end
            # @info "$(Threads.threadid())th thread : $(epoch)th epoch completed in $((Dates.now() - t_ini).value/1000) sec"
        end
    else
        for interval in 1:n_interval
            for i in 1:debug_interval
                epoch += 1
                for (x, y) in train_loader
                    gs = Flux.gradient(() -> loss(x, y), ps)
                    Flux.Optimise.update!(opt, ps, gs)
                    tip.matrix .= min.(tip.matrix, 0.0)
                    tip.matrix .= MDT.translate_tip_mean(tip.matrix)
                end
                loss_train[epoch] = loss(image_data_copy, image_data)
            end
            # @info "$(Threads.threadid())th thread : $(epoch)th epoch completed in $((Dates.now() - t_ini).value/1000) sec"
        end
    end

    if !need_loss_minimizing
        min_loss_tip = Surface(Matrix{Float32}(undef, 0, 0), mold=images[1])
    end

    DifferentiableBTRResult(lambda, max_epoch, tip, min_loss_tip, loss_train)
end

end