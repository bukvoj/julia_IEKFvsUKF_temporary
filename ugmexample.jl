N = 10000 # Number of simulations
Tmax = 50 # Number of time steps

σ_ω = 0.1 * ones(1,1)# Process noise
σ_v = 0.1 * ones(1,1) # Measurement noise

steplength = 0.5

outputpath = "ugm_benchmark_highnoise.png"

using LowLevelParticleFilters
using LinearAlgebra
using Random, Distributions
using Plots

Random.seed!(42) # Setting the seed

function rms(data::Vector{Float64})
    return sqrt(mean(data .^ 2))
end

include("models.jl")


GTx = zeros(Tmax,N)
UKFx = zeros(Tmax,N)
EKFx = zeros(Tmax,N)
IEKFx = zeros(Tmax,N)

dynamics(x,u,p,t) = ugmd(x,t)
measurement(x,u,p,t) = ugmm(x,t)

measnoise = Normal(0.0, σ_v[1])
procnoise = Normal(0.0, σ_ω[1])


for i in 1:N
    global xhat = 0.0
    global P = 1.0
    global ukf = UnscentedKalmanFilter(dynamics, measurement, σ_ω, σ_v, MvNormal([xhat],[P]); nu=1, ny=1, p=nothing)
    global ekf = ExtendedKalmanFilter(dynamics, measurement, σ_ω, σ_v, MvNormal([xhat],[P]); nu=1, ny=1, p=nothing)
    global iekf = IteratedExtendedKalmanFilter(dynamics, measurement, σ_ω, σ_v, MvNormal([xhat],[P]); nu=1, ny=1, p=nothing, step=steplength)

    x = 0.1

    for t in 1:Tmax

        x = ugmd(x, t) .+ rand(procnoise)
        global y = [ugmm(x, t) .+ rand(measnoise)]

        predict!(ukf,nothing,nothing,t)
        correct!(ukf, nothing,y,nothing,t)

        predict!(ekf,nothing,nothing,t)
        correct!(ekf, nothing,y,nothing,t)

        predict!(iekf,nothing,nothing,t)
        correct!(iekf, nothing, y,nothing,t)

        # save the error
        UKFx[t,i] = norm(ukf.x .- x)
        EKFx[t,i] = norm(ekf.x .- x)
        IEKFx[t,i] = norm(iekf.x .- x)
    end
end

ukferr = GTx - UKFx
ekferr = GTx - EKFx
iekferr = GTx - IEKFx
for t in 1:Tmax
    println("t = $t, RMSE = ", rms(ukferr[t,:]), " (UKF), ", rms(ekferr[t,:]), " (EKF)", ", ", rms(iekferr[t,:]), " (IEKF)")
end

ukfrms = [rms(ukferr[t,:]) for t in 1:Tmax]
ekfrms = [rms(ekferr[t,:]) for t in 1:Tmax]
iekfrms = [rms(iekferr[t,:]) for t in 1:Tmax]


plot(1:Tmax, ukfrms, label="UKF", xlabel="Time", ylabel="RMSE", title="RMSE vs Time", lw=2)
plot!(1:Tmax, ekfrms, label="EKF", lw=2)
plot!(1:Tmax, iekfrms, label="IEKF", lw=2)
savefig(outputpath)



using BenchmarkTools
@benchmark correct!(ukf, nothing,[y],nothing,10)
@benchmark correct!(iekf, nothing,[y],nothing,10)

@btime correct!(ukf, nothing,[y],nothing,10)
@btime correct!(iekf, nothing,[y],nothing,10)