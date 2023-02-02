module SPM

using Unitful
using CUDA

include("spm_core.jl")
include("hdr.jl")
include("btr.jl")
include("plots.jl")

using .SPMCore  # SPMをincludeすればSPMCore下

end
