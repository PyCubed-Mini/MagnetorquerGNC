function hat(ω::Vector)
    return [0 ω[3] -ω[2]
        -ω[3] 0 ω[1]
        ω[2] -ω[1] 0]
end

function randomMatrix(covariance)
    ϕ = √covariance * randn(3)
    return exp(hat(ϕ))
end