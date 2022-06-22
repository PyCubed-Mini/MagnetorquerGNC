# [Simulator/simulator.jl]

""" Simulator 

    Contains code for simulating a satellite. Includes position and orientation terms. 
"""
module Simulator 
    
    using SatelliteDynamics, LinearAlgebra 
    include("mag_field.jl")
    include("quaternions.jl")
    include("kf.jl")


    export initialize_orbit, intialize_orbit_oe   # Generate valid initial conditions for a satellite 
    export rk4      # Interface function for updating state 
    export IGRF13   # Function call that provides magnetic field vector given position and time
    export my_sim
    export my_sim_kf

    
    """
        dynamics(x, J, u, t; Rₑ, σβ)
        
          Propagates the state dynamics for a satellite, where the state is (probably) defined 
        as [position, velocity, scalar-first quaternion, angular velocity]
            'x = [r⃗, v⃗, (q₀, q⃗), ω]'

        Includes gravity from spherical harmonics, uniform drag, SRP, and third body 
        gravity from the sun and moon. 

        Arguments:
        - x:   STATE struct                                             |  [13,]
        - J:   Inertia matrix                                           |  [3, 3] 
        - u:   control input (for rotation only)                        |  [3,]    
        - t:   Current time (as an epoch)                               |  Epoch
        - Rₑ:  (Optional) Radius of the Earth (default is meters)       |  Scalar

        Returns:
        - x:   updated STATE struct                                     |  [13,]
    """
    function dynamics(x, J, u, t; Rₑ = 6378136.3, kwargs...)

        # Split apart the state (NOTE that this assumes you don't change any ordering...)
        r, v, q, ω = x[1:3], x[4:6], x[7:10], x[11:13]

        if norm(r) < Rₑ
            error("Error: Satellite impacted Earth at time $t!")
        end


        ṙ = v 
        v̇ = accel_perturbations(t, r, v; kwargs...)
        q̇ = qdot(q, ω)
        ω̇ = J \ (u - cross(ω, J * ω))
        
        return [ṙ; v̇; q̇; ω̇]
    end

    """
        accel_perturbations(epc, r, v; mass, area_drag, coef_drag, area_srp, coef_srp, n_grav, m_grav, third_body)

        Generates the acceleration for a spacecraft in LEO. Accounts for a variety of factors, 
        including spherical harmonics, atmospheric drag, SRP, and thirdbody from sun and moon

        ForwardDiff friendly (written by Kevin)
    """
    function accel_perturbations(epc::Epoch, r, v ;
                                mass::Real=1.0, area_drag::Real=0.01, coef_drag::Real=2.3, area_srp::Real=1.0, 
                                coef_srp::Real=1.8, n_grav::Integer=10, m_grav::Integer=10, third_body::Bool = true)

        # These functions don't like static, so we gotta adjust
        x = [r; v]

        # Compute ECI to ECEF Transformation -> IAU2010 Theory
        PN = bias_precession_nutation(epc)
        E  = earth_rotation(epc)
        W  = polar_motion(epc)
        R  = W * E * PN

        # Compute sun and moon position
        r_sun  = sun_position(epc)
        r_moon = moon_position(epc)

        # Compute acceleration (eltype(x) makes this forward diff friendly)
        a = zeros(eltype(x), 3)

        # spherical harmonic gravity
        a += accel_gravity(x, R, n_grav, m_grav)

        # atmospheric drag
        ρ = density_harris_priester(epc,r)
        a += accel_drag([r;v],ρ,mass, area_drag, coef_drag, Array{Real, 2}(PN))

        # SRP
        nu = eclipse_cylindrical(x, r_sun)  # conical doesn't work correctly
        a += nu * accel_srp(x, r_sun, mass, area_srp, coef_srp)

        if third_body
            # third body sun
            a += accel_thirdbody_sun(x, r_sun)

            # third body moon
            a += accel_thirdbody_moon(x, r_moon)
        end

        return (a)
    end

    """
        rk4(x, J, u, t, h)

        Modified RK4 function for integrating state. Normalizes quaternions before returning.
    """
    function rk4(x, J, u, t, dt; kwargs...)
        k₁ = dt * dynamics(x,        J, u, t       ; kwargs...)
        k₂ = dt * dynamics(x + k₁/2, J, u, t + dt/2; kwargs...)
        k₃ = dt * dynamics(x + k₂/2, J, u, t + dt/2; kwargs...)
        k₄ = dt * dynamics(x + k₃,   J, u, t + dt  ; kwargs...)

        # Normalize the quaternion because we aren't doing any fancy se(3) stuff 
        x⁺ = x + (1/6) * (k₁ + 2 * k₂ + 2 * k₃ + k₄)
        x⁺[7:10] /= norm(x⁺[7:10])

        return x⁺
    end

    """ initialize_orbit(; r, v, a, q, ω, Rₑ, μ)

          Initializes a random, viable orbit given a few different terms, usually 
        a position 'r' in Cartesian coordinates. Initial velocity may be specified, but 
        if specified it will not necessarily result in a stable orbit. 

        The initial starting position, velocity, semi-major axis, orientation, and angular 
        velocity may be either specified or determined randomly. 

        Arguments:
         - r:  (Optional) Height above ground that the satellite will start its orbit at    |  Scalar 
         - v:  (Optional) Magnitude of initial velocity                                     |  Scalar 
         - a:  (Optional) Semi-major axis                                                   |  Scalar 
         - q:  (Optional) Initial attitude, as a unit quaternion                            |  [4,]
         - ω:  (Optional) Initial angular velocity                                          |  [3,]

        Returns:
         - x:  Initial state, as 
                    x = [r, v, q, ω]
    """
    function initialize_orbit(; r = nothing, v = nothing, a = nothing, q = nothing, ω = nothing, Rₑ = 6378.1363e3, μ = 3.9860044188e14)

        ### POSITION ###

        # If unspecified, generate a random altitude for the satellite 
        _r = !isnothing(r) ? r : (rand(100:1000) * 1000) + Rₑ # Random height, in km 

        # If unspecified, pick some amount of eccentricity 
        _a = !isnothing(a) ? a * _r : rand(1.0:0.05:1.5) * _r   # Semi-major axis. Results in a circular orbit for a == r

        # If unspecified, calculate the necessary initial velocity to maintain a valid orbit 
        _v = !isnothing(v) ? v : sqrt(μ*( (2/_r) - (1/_a)))

        # Now we have to pick directions for r and v... 
        axes = [ [1, 0, 0], [0, 1, 0], [0, 0, 1] ]
        ax   = rand(1:3)
        r₀ = _r * splice!(axes, ax)  # removes the selected axis from axes 
        v₀ = _v * rand(axes)         # Remaining options are both perpendicular and valid 

        #   allow for negative axes...
        r₀ = (randn() > 0.0)  ?  r₀ : -r₀ 
        v₀ = (randn() > 0.0)  ?  v₀ : -v₀


        ### ORIENTATION ###

        # If unspecified, generate a random starting orientation
        q₀ = !isnothing(q) ? q : randn(4);
        q₀ /= norm(q₀)  # Ensure it is unit 

        # If unspecified, generate a random starting angular velocity
        ω₀ = !isnothing(ω) ? ω : 0.5 * randn(3)

        x = [r₀; v₀; q₀; ω₀]
        return x
    end

    """ initialize_orbit_oe(; ϵ, i, Ω, ω, M, r, Rₑ, detumbled, use_degrees)

          Initializes a random, viable orbit given a set of orbital elements. 
        These orbital elements can be specified, but if unspecified they will be somewhat 
        similar to the ISS orbit (I think...) with some amount of random perturbations each time.

        Arguments:
         - [ϵ, i, Ω, ω, M]:  (Optional) Orbital elements                                    |  Scalars
         - r:    (Optional) Height of satellite (above ground, in m)                        |  Scalar 
         - Rₑ:   (Optional) Radius of the Earth (in m)                                      |  Scalar 
         - detumbled:  (Optional) Flat that adjusts how much the satellite is spinning      |  Bool
                            (NOTE that detumbled just means 𝑙𝑜𝑤 angular velocity, not zero)

        Returns: 
         - x:   Vector containing inital state                                              |  [13,]
                    x = [r, v, q, ω]
         - T_orbit:  Amount of seconds required for a full orbit given x                    |  Scalar
    """
    function intialize_orbit_oe(; ϵ = 0.0001717 + 0.00001 * randn(), i = 51.6426 + randn(), Ω = 178.1369 + randn(), 
                                  ω = 178.1369 + randn(), M = 330.7918 + 100 * randn(), r = 421e3, Rₑ = 6378136.3, 
                                  detumbled = false) 
        # ecc = 0.0001717 + 0.00001 *randn()
        # inc = 51.6426 + randn()
        # Ω   = 178.1369 + randn()
        # ω   = 174.7410 + randn()
        # M   = 330.7918 + 100 + randn()   # +94/95 is just before sun, -40 is just before eclipse
        sma = (Rₑ + r) / (1 + ecc)  # Apogee = semi_major * (1 + ecc)
    
        oe0 = [sma, ecc, inc, Ω, ω, M]   # Initial state, oscullating elements
        eci0 = sOSCtoCART(oe0, use_degrees = true) # Convert to Cartesean
    
        r₀ = (eci0[1:3])
        v₀ = (eci0[4:6])
        q₀ = randn(4);  q₀ = SVector{4, Float64}(q₀ / norm(q₀))
        ω₀ = (detumbled) ? (0.05 * randn(3)) : (0.5 * randn(3))
        
        T_orbit = orbit_period(oe0[1])
        x = [r₀; v₀; q₀; ω₀]
    
        return x, T_orbit 
    end

    function my_sim(control_fn)
        x₀ = initialize_orbit() 
        println("intialized orbit!")
        # x₀[11:13] .=0
        # x₀[11:13] *= x₀[11:13].*10.0 # Spinning very fast

        J  = [0.3 0 0; 0 0.3 0; 0 0 0.3]  # Arbitrary inertia matrix for the Satellite 
        t  = Epoch(2020, 11, 30)          # Starting time is Nov 30, 2020
        dt = 0.5                          # Time step, in seconds

        N = 1000000
        q_hist = zeros(N, 4)
        q_hist[1, 1:3] .= x₀[11:13]
        q_hist[1, 4] = norm(x₀[11:13])
        x = x₀
        for i = 1:N - 1
            r, v, q, ω = x[1:3], x[4:6], x[7:10], x[11:13]
            b = IGRF13(r, t)
            x = rk4(x, J, control_fn(ω, b), t, dt)
            t += dt                      # Don't forget to update time (not that it really matters...)
            # q_hist[i + 1, :] .= x[7:10]
            q_hist[i + 1, 1:3] .= x[11:13]
            ω = norm(x[11:13])
            if ω < 0.001
                q_hist = q_hist[1:i, :]
                break
            end
            q_hist[i + 1, 4] = ω
        end

        return q_hist

    end

    function my_sim_kf(control_fn)
        x₀ = initialize_orbit() 
        println("intialized orbit!")
        # x₀[11:13] .=0
        # x₀[11:13] *= x₀[11:13].*10.0 # Spinning very fast

        J  = [0.3 0 0; 0 0.3 0; 0 0 0.3]  # Arbitrary inertia matrix for the Satellite 
        t  = Epoch(2020, 11, 30)          # Starting time is Nov 30, 2020
        dt = 0.5                          # Time step, in seconds

        N = 1000000
        q_hist = zeros(N, 9)
        q_hist[1, 1] = norm(x₀[11:13])
        x = x₀
        
        q₀ = x₀[7:10]
        P = Diagonal([1, 1, 1, 1])
        kf = EKF(q₀, P)

        for i = 1:N - 1
            r, v, q, ω = x[1:3], x[4:6], x[7:10], x[11:13]
            b = IGRF13(r, t)
            x = rk4(x, J, control_fn(ω, b), t, dt)
            t += dt                      # Don't forget to update time (not that it really matters...)
            q_hist[i + 1, 1] = norm(ω)
            q_hist[i + 1, 2:5] .= q

            q_hist[i+1, 6:9] .= step(kf, q)

            if ω < 0.001
                q_hist = q_hist[1:i, :]
                break
            end
        end

        return q_hist

    end
end