using Quaternions

struct EKF
    q::Quaternion
    P::Matrix{Float64}
end

function step(e::EKF, qn::Tuple{Float64,Float64,Float64})
    e.q = qn
    return e.q
end