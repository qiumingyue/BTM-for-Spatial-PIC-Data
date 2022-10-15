using Distributions, LinearAlgebra, Random, NLsolve, CSV, DataFrames, StatsPlots, LaTeXStrings, MCMCChains

cd("E:\\Spatial\\github-spatial")
include("Functions.jl")

W = Matrix(CSV.read("adjacency matrix.csv",DataFrame,header=1))
data = Matrix(CSV.read("simulated data--PH.csv",DataFrame))
r = 0
n_knots = 15
order_1 = 3
a_eta = b_eta = 1
a_tau = b_tau = 1
niter = 14000
global sigma = 10
global MH_iter = 14000
global c0 = 1.5
global MH_iter_spatial = 14000
global c1 = 1.5
global g = 1
est(r, order_1, n_knots, niter, a_eta, b_eta, a_tau, b_tau, data)

