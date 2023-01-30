module SPM

using Unitful
using CUDA

include("spm_core.jl")
include("hdr.jl")
include("btr.jl")
include("plot.jl")

using .SPMCore  # SPMをincludeすればSPMCore下

end
