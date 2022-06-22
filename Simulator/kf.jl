# [Simulator/kf.jl]

using Quaternions

struct EKF
    q::Vector{Float64}
    P::Matrix{Float64}
end

function step(e::EKF, qn::Vector{Float64})
    e.q = qn
    return e.q
end