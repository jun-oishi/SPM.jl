
cd(@__DIR__)
import Pkg
Pkg.activate("../")
print("loading SPM.jl...")
include("../src/SPM.jl")
println("done")

src = "B3446"
print("loading data/$(src).hdr...")
image = SPM.HDR.loadHDR("data/$(src).hdr")
println("done")

extracted = Vector{SPM.Surface}(undef, 64)
for i=0:7, j=0:7
    extracted[i*8+j+1] = SPM.extract(image, (i*64+1, j*64+1), (64, 64))
end

print("running BTR...")
result = SPM.BTR.solveDifferentiableBTR(extracted, 10, 100, 0.1)
println("done")

SPM.BTR.saveResult(result, "data/$(src)_btr_result")

refined = SPM.BTR.refineImage(image, result.tip)

SPM.HDR.saveHDR(refined, "data/$(src)_refined.hdr")
println("results saved")