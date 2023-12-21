module SPM

include("spm_core.jl")
include("hdr.jl")
# include("btr.jl")
include("plots.jl")

using .SPMCore  # SPMをincludeすればSPMCore下で定義されたものはSPM.で呼び出せる
Plots = SPMPlots # conflictを避けるためにSPMPlotsで定義したが長いのでPlotsにalias

end
