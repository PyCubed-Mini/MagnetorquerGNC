using Random, Distributions

mutable struct Gyro
    β::Vector{Float64}
    Vβ::Matrix{Float64}
    Vω::Matrix{Float64}
end

""" update(g, ω)
    Arguments:
        - g: The gyro                     | Gyro
        - ω: True angular velocity        | [3,]
    Returns:
        - ω: Measured angular velocity    | [3,]

"""
function update(g::Gyro, ω::Vector{Float64})
    g.β += rand(MvNormal([0; 0; 0], g.Vβ))
    return ω + g.β + rand(MvNormal([0; 0; 0], g.Vω))
end