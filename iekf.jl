using ForwardDiff
using LinearAlgebra

# Implements the data (measurement, correction) step of the iekf algorithm

# This is to differentiate between scalar and multidimensional system and measurement
for XTYPE in ((Vector{<:Real},AbstractArray{<:Real}),(Real,Real))
for YTYPE in ((Vector{<:Real},AbstractArray{<:Real}),(Real,Real))
@eval begin

function iekfcorrect(h::Function, 
                    x::$XTYPE[1],
                    P::$XTYPE[2],
                    y::$YTYPE[1],
                    R::$YTYPE[2],
                    maxiters=30,eps=1e-8,steplen=1.0,
                    u...)
    xi = Float64.(x)

    # function V(x,y,R,P,xhat, h)
    #     (y-h(x,u...))'/R*(y-h(x,u...)) + (xhat-x)'/P*(xhat-x)
    # end

    hh = args -> h(args,u...)
    Hi = zeros(length(y),length(x))

    i = 1
    y_hat = similar(y)
    while true
        prev = xi
        ForwardDiff.jacobian!(Hi, hh, xi)
        y_hat = h(xi,u...)
        residual =  y - y_hat
        Si = Hi*P*Hi' + R
        Ki = P*Hi' *(Si\I) 

        # VARIABLE STEP LENGTH
        Δ = x - xi + Ki*(residual-Hi*(x-xi))
        xi = x + steplen*Δ
        # # FULL STEP
        # xi = x + Ki*(residual-Hi*(x-xi))
        if sum(xi-prev) < eps || i == maxiters
            return xi, P-Ki*Hi*P, residual, Si, Ki, Hi, i
        end
        i += 1
    end
end


end #eval
end #for
end #for



# This file contains the predict step of the EKF and IEKF algorithms.

# This is to differentiate between scalar and multidimensional system
for T in ((Vector{<:Real},AbstractArray{<:Real}),(Real,Real)) 
@eval begin

function iekfpredict(f::Function, 
                    x::Union{Vector, Real},
                    P::$T[2],
                    Q::$T[2],
                    u...)
    ff = args -> f(args,u...)
    F = ForwardDiff.jacobian(ff, Float64.(x))
    x_new = f(x,u...)
    P_new = F[1]*P*F[1]' + Q
    return x_new,P_new
end

end #eval
end #for