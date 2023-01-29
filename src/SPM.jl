module SPM

# Write your package code here.

include("spm_core.jl")
include("btr.jl")
include("hdr.jl")
include("plot.jl")

using .SPMCore  # SPMをincludeすればSPMCore下

end
