module SPM

using Unitful

include("spm_core.jl")
include("hdr.jl")
include("btr.jl")
include("plot.jl")

using .SPMCore  # SPMをincludeすればSPMCore下

end
