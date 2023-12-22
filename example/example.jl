
import Pkg
Pkg.activate("../")
import Plots
include("../src/SPM.jl")

# データの読み込み
image = SPM.loadImage("data/")