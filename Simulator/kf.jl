# [Simulator/kf.jl]

using LinearAlgebra
include("quaternions.jl")

mutable struct EKF
    q::Vector{Float64}
    β::Vector{Float64}
    P::Matrix{Float64}
end

function hat(ω::Vector)
    return [0 ω[3] -ω[2]
        -ω[3] 0 ω[1]
        ω[2] -ω[1] 0]
end

function f(
    q::Vector{Float64},
    β::Vector{Float64},
    ω::Vector{Float64},
    δt::Float64
)
    θ = norm(ω - β) * δt
    r = (ω - β) / norm(ω - β)
    return L(q) * [cos(θ / 2); r * sin(θ / 2)]
end

function step(
    e::EKF,
    ω::Vector{Float64},
    δt::Float64,
    ⁿr_mag::Vector{Float64},
    ⁿr_sun::Vector{Float64},
    ᵇr_mag::Vector{Float64},
    ᵇr_sun::Vector{Float64}
)
    # Predict:
    qₚ = f(e.q, e.β, ω, δt) # β remains constant
    R = exp(-hat(ω - e.β) * δt)
    A = [
        R (-δt*I(3))
        zeros(3, 3) I(3)
    ]
    W = I(6) * 0.01 # related to some noise or something (ask Zac)
    Pₚ = A * e.P * A' + W

    # Innovation
    Q = quaternionToMatrix(qₚ)
    Z = [ᵇr_mag; ᵇr_sun] - [Q zeros(3, 3); zeros(3, 3) Q] * [ⁿr_mag; ⁿr_sun]
    C = [hat(ᵇr_mag) zeros(3, 3); hat(ᵇr_sun) zeros(3, 3)]
    V = I(6) * 0.01 # Something else
    S = C * Pₚ * C' + V

    # Kalman Gain
    L = Pₚ * C' * inv(S)

    # Update
    println("Z: ", Z)
    δx = L * Z
    println(size(L), size(Z))
    println(δx, size(δx))
    ϕ = δx[1:3]
    δβ = δx[4:6]
    θ = norm(ϕ)
    r = ϕ / θ
    qᵤ = ⊙(qₚ, [cos(θ / 2); r * sin(θ / 2)])
    βᵤ = e.β + δβ
    Pᵤ = (I(6) - L * C) * Pₚ * (I(6) - L * C)' + L * V * L'

    e.q = qᵤ
    e.β = βᵤ
    e.P = Pᵤ
    println("β: ", e.β, "δβ: ", δβ)
    return e
end