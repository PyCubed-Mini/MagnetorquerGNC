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
    r = (ω - β) / norm(ω - e.β)
    return [L(q) * [r .* sin(θ / 2); cos(θ / 2)], β]
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
    xₚ = f(e.q, e.β, ω, δt)
    R = exp(-hat(ω - e.β) * δt)
    A = [
        R (-δt*I(3))
        zeros(3, 3) I(3)
    ]
    W = zeros(6, 6) # related to some noise or something (ask Zac)
    Pₚ = A * P * A' + W
    # Innovation
    Z = [ᵇr_mag; ᵇr_sun] - [Q zeros(3, 3); zeros(3, 3) Q] * [ⁿr_mag; ⁿr_sun]
    C = [hat(ᵇr_mag) zeroes(3, 3); hat(ᵇr_sun) zeros(3, 3)]
    V = zeros(3, 3) # Something else
    S = C * Pₚ * C' + V
    # Kalman Gain
    L = Pₚ * C' * inv(S)
    # Update
    δx = L * Z
    return e.q
end