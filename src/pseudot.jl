module pseudot
    using KernelFunctions: Kernel
    # using KernelFunctions
    # using Mamba
    export Trejactory, Trace, bifurcateKernel
    struct Trejactory
        t::Float64
        c::Int64
        xr::Float64
        Trejactory(t, c, xr) = !(0<=t<=1) | !(1<=c<=2) | !(-1<xr<1) ? error("invalid input") : new(t, c, xr)
    end
    dt(x::Trejactory, t::Float64) = (x.t - t)
    dt2(x::Trejactory, t::Float64) = (x.t - t)^2
    dT(x::Trejactory, y::Trejactory) = (x.t - y.t)
    dT2(x::Trejactory, y::Trejactory) = (x.t - y.t)^2
    dx2(x::Trejactory, y::Trejactory) = (x.xr - y.xr)^2
    Txt(x::Trejactory, t::Float64) = x.t * t
    TxT(x::Trejactory, y::Trejactory) = x.t * y.t
    TxTr(x::Trejactory, y::Trejactory) = x.xr * y.xr
    Tadj(x::Trejactory, t::Float64) = Trejactory(dt(x, t), x.c, x.xr)
    Phase1(x1::Trejactory, x2::Trejactory, t::Float64) = (x1.t < t) & (x2.t < t)
    Phase2(x1::Trejactory, x2::Trejactory, t::Float64) = (x1.t > t) & (x2.t > t)
    SameType(x1::Trejactory, x2::Trejactory) = x1.c == x2.c
    struct Trace
        X::Vector{Trejactory}
    end
    Trace(; X::Vector{Trejactory}) = Trace(X)
    Base.length(L::Trace) = length(L.X)
    Tt(L::Trace) = [x.t for x in L.X]
    root(L::Trace, Tb::Float64) = [x.t <= Tb for x in L.X]
    branch1(L::Trace, Tb::Float64) = [(x.t > Tb) & (x.c == 1) for x in L.X]
    branch2(L::Trace, Tb::Float64) = [(x.t > Tb) & (x.c == 2) for x in L.X]
    struct bifurcateKernel <: Kernel
        Tb::Float64
        lambda::Vector{Float64}  
        sigma::Float64
        a::Float64 
        b::Int64
    end
    bifurcateKernel(; Tb::Float64, lambda::Vector{Float64}, sigma::Float64, a::Float64, b::Int64) = bifurcateKernel(Tb, lambda, sigma, a, b) 
    function (k::bifurcateKernel)(x::Trejactory, y::Trejactory)
        if SameType(x, y)
            if Phase1(x, y, k.Tb)
                K = k.lambda[1] * exp(-k.sigma * dT2(x, y)) + (k.lambda[2] * TxT(x, y) + k.a) ^ k.b 
            elseif Phase2(x, y, k.Tb)
                K = k.lambda[1] * exp(-k.sigma * dT2(x, y)) + (k.lambda[2] * TxT(Tadj(x, k.Tb), Tadj(y, k.Tb)) + k.a) ^ k.b 
            else 
                K = k.lambda[1] * exp(-k.sigma * dt2(x, k.Tb)) * exp(-k.sigma * dt2(y, k.Tb)) 
            end
        else
            if Phase1(x, y, k.Tb)
                K = k.lambda[1] * exp(-k.sigma * dT2(x, y)) + (k.lambda[2] * TxT(x, y) + k.a) ^ k.b
            else 
                K = k.lambda[1] * exp(-k.sigma * dt2(x, k.Tb)) * exp(-k.sigma * dt2(y, k.Tb))
            end
        end
        #return  K + (k.lambda[3] * TxTr(x, y) + k.a) ^ k.b
        return K + k.lambda[3] * exp(-1e5 * dx2(x, y))
    end

end


