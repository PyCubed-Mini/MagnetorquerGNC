function randomMatrix(covariance)
    ϕ = √covariance * randn(3)
    return exp(hat(ϕ))
end