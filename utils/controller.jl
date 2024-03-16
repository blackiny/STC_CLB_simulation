using LinearAlgebra, Plots
import ForwardDiff as FD
using SparseArrays
import MeshCat as mc 
using Test
using Random
import Convex as cvx 
import ECOS
using ProgressMeter

"""
In each function with state x as an argument, x refers to deviation from equilibrium points
"""
function V(params::NamedTuple, s::Float64)
    v_max, s_st, s_go = params.v_max, params.s_st, params.s_go
    v = 0.0
    if s > s_st && s ≤ s_go
        v = 0.5 * v_max * (1 - cos(π*((s-s_st)/(s_go-s_st))))
    elseif s ≥ s_go
        v = v_max
    end
    return v
end

function V_devirative(params::NamedTuple, s::Float64)
    v_max, s_st, s_go = params.v_max, params.s_st, params.s_go
    ∂v = 0.0
    if s > s_st && s ≤ s_go
        ∂v = 0.5 * v_max * sin(π*((s-s_st)/(s_go-s_st))) * π / (s_go-s_st)
    end
    return ∂v
end

function r_steady(params::NamedTuple, t::Float64, speed::Float64)
    return speed
end

function r_brake_acc(params::NamedTuple, t::Float64, speed::Float64, acc::Float64, t_brake::Float64, t_acc::Float64)
    s = speed
    if t > t_brake && t ≤ t_acc
        s = max(0.0, speed - acc*(t-t_brake)) 
    elseif t > t_acc && t ≤ 2*t_acc - t_brake
        s_min = max(0.0, speed - acc*(t_acc-t_brake))
        s = min(speed, s_min + acc*(t-t_acc)) 
    end
    return s
end

function naive_observer(x::Vector)
    return x
end

function string_linear_dynamics(params::NamedTuple, x::Vector, u::Float64, r, t::Float64)
    # u is zero-order Hold
    # r is a function of t
    # f, g meet the eq: ẋ = f + g*u 
    N = length(x) ÷ 2 - 1 #index starting from 0
    a, b = params.a, params.b
    s_star = params.s_star
    a1, a2, a3 = a*V_devirative(params, s_star), a+b, b
    P0 = sparse([0. -1;0 0])
    P1 = sparse([0. -1;a1 -a2])
    Q1 = sparse([0. 1;0 a3])
    b0 = [0.;1]
    d0 = [1.;0]
    A = blockdiag(P0, kron(I(N), P1))
    for i=2:(N+1)
        A[2*(i-1).+(1:2), 2*(i-2).+(1:2)] .= Q1
    end
    B = [b0;zeros(2*N)]
    D = [d0;zeros(2*N)]
    f = A*x + D*r(params, t)
    g .= B
    return f + g*u, f, g
end

function string_linear_dynamics_(params::NamedTuple, x::Vector, u::Float64, r, t::Float64)
    ẋ, _, _ = string_linear_dynamics(params, x, u, r, t)
    return ẋ
end

function HDV_v_derivative(params::NamedTuple, x::Vector)
    N = length(x) ÷ 2 - 1
    a, b = params.a, params.b
    s_star, v_star = params.s_star, params.v_star
    s = x[1:2:end] .+ s_star
    v = x[2:2:end] .+ v_star
    ṡ = -diff(v)
    v̇ = zeros(N)
    for i=1:N
        v̇[i] = a*(V(params, s[i+1])-v[i+1]) + b*ṡ[i]
    end
    return v̇
end

function string_dynamics(params::NamedTuple, x::Vector, u::Float64, r, t::Float64)
    # u is zero-order Hold
    # r is a function of t
    # f, g meet the eq: ẋ = f + g*u 
    N = length(x) ÷ 2 - 1 #index starting from 0
    a, b = params.a, params.b
    f, g = zeros(length(x)), zeros(length(x))
    g[2] = 1.0
    f[1] = r(params, t) - x[2]
    f[3:2:end] .= -diff(x[2:2:end])
    f[4:2:end] .= HDV_v_derivative(params, x)
    return f + g*u, f, g
end

function string_dynamics_(params::NamedTuple, x::Vector, u::Float64, r, t::Float64)
    ẋ, _, _ = string_dynamics(params, x, u, r, t)
    return ẋ
end

function rk4(dynamics, params::NamedTuple, x::Vector, u::Float64, r, t0::Float64, dt::Float64)
    # vanilla RK4
    k1 = dt*dynamics(params, x, u, r, t0)
    k2 = dt*dynamics(params, x + k1/2, u, r, t0+dt/2)
    k3 = dt*dynamics(params, x + k2/2, u, r, t0+dt/2)
    k4 = dt*dynamics(params, x + k3, u, r, t0+dt)
    x + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
end

function nominal_controller(params::NamedTuple, x::Vector, r, t::Float64)
    N = length(x) ÷ 2 - 1 #index starting from 0
    a, b = params.a, params.b
    s_star = params.s_star
    a1, a2, a3 = a*V_devirative(params, s_star), a+b, b
    μ, k = params.μ, params.k
    @assert length(μ) == length(k) && length(k) == N
    u = a1*x[1]-a2*x[2]+a3*r(params, t)
    for i=1:N
        u += μ[i]*x[2*i+1] + k[i]*x[2*i+2]
    end
    return u
end

function CLB_TH(params::NamedTuple, x::Vector, r, t::Float64)
    N = length(x) ÷ 2 - 1
    s_star, v_star = params.s_star, params.v_star
    τ = params.τ
    s, v = x[1:2:end] .+ s_star, x[2:2:end] .+ v_star
    h = s - τ*v
    return h
end

function CLB_TH_derivative(params::NamedTuple, x::Vector, r, t::Float64)
    N = length(x) ÷ 2 - 1
    τ = params.τ
    ∂h = zeros(N+1, length(x))
    ∂h[1, 1] = 1.0
    ∂h[1, 2] = -τ
    for i=1:N
        ∂h[i+1, 1] = -1.0
        ∂h[i+1, 2] = τ
        ∂h[i+1, 2*i+1] = 1.0
        ∂h[i+1, 2*i+2] = -τ
    end
    return ∂h
end

function CLB_SDH(params::NamedTuple, x::Vector, r, t::Float64)
    N = length(x) ÷ 2 - 1
    s_star, v_star = params.s_star, params.v_star
    a = abs(params.u_min)
    τ = params.τ
    h = zeros(N+1)
    s, v = x[1:2:end] .+ s_star, x[2:2:end]
    h[1] = s[1] - τ*(v[1]-r(params, t)) - (v[1]-r(params, t))^2 / (2.0*a)
    v_diff = diff(v)
    h[2:end] .= s[2:end] - τ*v_diff - v_diff.^2 / (2.0*a)
    return h
end

function CLB_SDH_derivative(params::NamedTuple, x::Vector, r, t::Float64)
    N = length(x) ÷ 2 - 1
    s_star, v_star = params.s_star, params.v_star
    a = abs(params.u_min)
    τ = params.τ
    s, v = x[1:2:end] .+ s_star, x[2:2:end]
    ∂h = zeros(N+1, length(x))
    ∂h[1, 1] = 1.0
    ∂h[1, 2] = -(τ + (v[1]-r(params, t))/a)
    for i=1:N
        ∂h[i+1, 1] = -1.0
        ∂h[i+1, 2] = (τ + (v[1]-r(params, t))/a)
        ∂h[i+1, 2*i] += (τ + (v[i+1]-v[i])/a)
        ∂h[i+1, 2*i+1] = 1.0
        ∂h[i+1, 2*i+2] = -(τ + (v[i+1]-v[i])/a)
    end
    return ∂h
end
        
    

function CLB_controller(dynamics, CLB, CLB_derivative, params::NamedTuple, x::Vector, u_nominal::Float64, r, t::Float64)
    N = length(x) ÷ 2 - 1
    U = cvx.Variable(1)
    σ = cvx.Variable(N)

    p, γ = params.p, params.γ
    obj = cvx.sumsquares(U-u_nominal) + p*cvx.sumsquares(σ)
    prob = cvx.minimize(obj)
    prob.constraints += (σ ≥ 0)
    _, f, g = dynamics(params, x, u_nominal, r, t)
    h = CLB(params, x, r, t)
    ∂h = CLB_derivative(params, x, r, t)
    prob.constraints += (dot(∂h[1,:], f)+dot(∂h[1,:], g)*U+γ*h[1] ≥ 0)
    for i=1:N
        prob.constraints += (dot(∂h[i+1,:], f)+dot(∂h[i+1,:], g)*U+γ*h[i+1]+σ[i] ≥ 0)
    end
    cvx.solve!(prob, ECOS.Optimizer; silent_solver = true)
    σ = σ.value
    U = U.value

    return U, σ
end