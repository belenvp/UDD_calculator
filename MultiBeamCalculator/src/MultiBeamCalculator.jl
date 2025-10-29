module MultiBeamCalculator

const newaxis = [CartesianIndex()]
const c_light= 299792458          # Light Speed m/s
const h_planck = 4.13566733e-15   # eV
const hc = 1.2398e-06
const r_e = 2.8179403267e-15
const rho = 5e-4

function squeeze(array)
    dims = tuple(findall(==(1), size(array))...)
    dropdims(array; dims)
end


include("elements.jl")

import FFTW
include("setup.jl")

import LinearAlgebra: mul!
import Statistics: mean
import CUDA: CuArray, CUFFT
import KernelAbstractions as KA
import KernelAbstractions: @kernel, @index
import ProgressMeter: Progress, next!
include("simulation.jl")

import StatsBase: fit, uweights, Histogram, Weights
include("scans.jl")

end # module MultiBeamCalculator
