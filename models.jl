function ugmd(x,t)
    0.5 .*x .+ 25 .*x ./ (1 .+ x.^2) .+ 8 .*cos(1.2.*(t-1))
end

function ugmm(x,t)
    x.^2 ./20
end

function botd(x,t)
    x
end

function botm(x,t)
    [atan((x[2]-1.5)/(x[1]-0)), atan((x[2]-0)/(x[1]-0))]
end