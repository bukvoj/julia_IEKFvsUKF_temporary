N = 10000 # Number of simulations
Tmax = 50 # Number of time steps

σ_ω = 0.1 # Process noise
σ_v = 0.1 # Measurement noise

steplength = 0.5

outputpath = "ugm_benchmark_highnoise.png"

using LowLevelParticleFilters
using IteratedExtendedKalmanFilter
using Random, Distributions
using Plots

Random.seed!(42) # Setting the seed
include("iekf.jl")

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

measnoise = Normal(0.0, σ_v)
procnoise = Normal(0.0, σ_ω)


for i in 1:N
    global xhat = 0.0
    global P = 1.0
    global ukf = UnscentedKalmanFilter(dynamics, measurement, σ_ω*ones(1,1), σ_v*ones(1,1), MvNormal([xhat],[P]); nu=1, ny=1, p=nothing)
    ekf = ExtendedKalmanFilter(dynamics, measurement, σ_ω*ones(1,1), σ_v*ones(1,1), MvNormal([xhat],[P]); nu=1, ny=1, p=nothing)

    x = 0.1

    global y = 0.0

    for t in 1:Tmax
        x = ugmd(x, t) + rand(procnoise)
        y = ugmm(x, t) + rand(measnoise)

        predict!(ukf,nothing,nothing,t)
        correct!(ukf, nothing,[y],nothing,t)

        predict!(ekf,nothing,nothing,t)
        correct!(ekf, nothing,[y],nothing,t)

        
        # SORRY FOR THIS, IEKF EXPECTS VECTORS AND MATRICES, THIS SHOULD BE FIXED IN THE FINAL IMPLEMENTATION
        xhat, P = iekfpredict(ugmd, [xhat], [P]', [σ_ω], t)
        xhat = xhat[1]
        P = P[1,1]
        xhat, P = iekfcorrect(ugmm, [xhat], Matrix([P]'), [y], Matrix([σ_v]'), 10,1e-8, steplength, t)
        xhat = xhat[1]
        P = P[1,1]


        # save the ground truth, UKF estimate, EKF estimate, and IEKF estimate
        GTx[t,i] = x
        UKFx[t,i] = ukf.x[1]
        EKFx[t,i] = ekf.x[1]
        IEKFx[t,i] = xhat
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
@benchmark iekfcorrect(ugmm, [xhat], Matrix([P]'), [y], Matrix([σ_v]'), 10,1e-8, steplength, 10)
@benchmark correct!(ukf, nothing,[y],nothing,10)
@btime iekfcorrect(ugmm, [xhat], Matrix([P]'), [y], Matrix([σ_v]'), 10,1e-8, steplength, 10)
@btime correct!(ukf, nothing,[y],nothing,10)