N = 10000 # Number of simulations
Tmax = 25 # Number of time steps

σ_ω = 0.1 * [1.0 0.0; 0.0 1.0] # Process noise Q
σ_v = 1e-4 * [1.0 0.0; 0.0 1.0]# Measurement noise R

steplength = 0.5

outputpath = "bot_benchmark_lownoise.png"


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

dynamics(x,u,p,t) = botd(x,t)
measurement(x,u,p,t) = botm(x,t)

procnoise = MvNormal([0.0, 0.0], σ_ω)
measnoise = MvNormal([0.0, 0.0], σ_v)

xhat = [1.5, 1.5]
P = 0.1 * [1.0 0.0; 0.0 1.0]
for i in 1:N
    global xhat = [1.5, 1.5]
    global P = 0.1 * [1.0 0.0; 0.0 1.0]
    global ukf = UnscentedKalmanFilter(dynamics, measurement, σ_ω, σ_v, MvNormal(xhat,P); nu=1, ny=2, p=nothing)
    global ekf = ExtendedKalmanFilter(dynamics, measurement, σ_ω, σ_v, MvNormal(xhat,P); nu=1, ny=2, p=nothing)
    global iekf = IteratedExtendedKalmanFilter(dynamics, measurement, σ_ω, σ_v, MvNormal(xhat,P); nu=1, ny=2, p=nothing, step=steplength)

    x = [1.5,1.5]

    for t in 1:Tmax
        x = botd(x, t) + rand(procnoise)
        global y = botm(x, t) + rand(measnoise)

        predict!(ukf,nothing,nothing,t)
        correct!(ukf, nothing,y,nothing,t)

        predict!(ekf,nothing,nothing,t)
        correct!(ekf, nothing,y,nothing,t)

        predict!(iekf,nothing,nothing,t)
        correct!(iekf, nothing, y,nothing,t)

        # save the error
        UKFx[t,i] = norm(ukf.x - x)
        EKFx[t,i] = norm(ekf.x - x)
        IEKFx[t,i] = norm(iekf.x - x)
    end
end

ukferr = UKFx
ekferr = EKFx
iekferr = IEKFx
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
@benchmark correct!(ukf, nothing,y,nothing,10)
@benchmark correct!(iekf, nothing,y,nothing,10)

@btime correct!(ukf, nothing,y,nothing,10)
@btime correct!(iekf, nothing,y,nothing,10)