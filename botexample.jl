N = 10000 # Number of simulations
Tmax = 25 # Number of time steps

σ_ω = 1 * [1.0 0.0; 0.0 1.0] # Process noise Q
σ_v = 1e-4 * [1.0 0.0; 0.0 1.0]# Measurement noise R

steplength = 0.5

outputpath = "bot_benchmark_highnoise.png"


using LowLevelParticleFilters
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
    ekf = ExtendedKalmanFilter(dynamics, measurement, σ_ω, σ_v, MvNormal(xhat,P); nu=1, ny=2, p=nothing)

    x = [1.5,1.5]

    for t in 1:Tmax
        x = botd(x, t) + rand(procnoise)
        global y = botm(x, t) + rand(measnoise)

        predict!(ukf,nothing,nothing,t)
        correct!(ukf, nothing,y,nothing,t)

        predict!(ekf,nothing,nothing,t)
        correct!(ekf, nothing,y,nothing,t)

        xhat, P = iekfpredict(botd, xhat, P, σ_ω, t)
        xhat, P = iekfcorrect(botm, xhat, P, y, σ_v, 10,1e-8,steplength, t)

        # save the error
        UKFx[t,i] = sqrt((ukf.x[1]-x[1])^2 + (ukf.x[2]-x[2])^2)
        EKFx[t,i] = sqrt((ekf.x[1]-x[1])^2 + (ekf.x[2]-x[2])^2)
        IEKFx[t,i] = sqrt((xhat[1]-x[1])^2 + (xhat[2]-x[2])^2)
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
@benchmark iekfcorrect(botm, xhat, P, y, σ_v, 10,1e-8,steplength, 10)
@benchmark correct!(ukf, nothing,y,nothing,10)
@btime iekfcorrect(botm, xhat, P, y, σ_v, 10,1e-8,steplength, 10)
@btime correct!(ukf, nothing,y,nothing,10)