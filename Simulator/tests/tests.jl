# [Simulator/test/tests.jl]

using Test, LinearAlgebra, SatelliteDynamics, Plots, JSON
include("../simulator.jl");   using .Simulator


@testset "General Function Call" begin 
    """ Test general function calls and ensure that everything looks as it should """
    # println("General Function Call")
    x₀ = initialize_orbit() 

    # Ensure that r and v are orthogonal 
    @test dot(x₀[1:3], x₀[4:6]) ≈ 0.0 
    
    # Ensure that the quaternion is unit 
    @test norm(x₀[7:10]) ≈ 1.0

    J  = [0.2 0 0; 0 0.2 0; 0 0 0.2]  # Arbitrary inertia matrix for the Satellite 
    t  = Epoch(2021, 5, 31)           # Starting time is May 31, 2021
    dt = 1.0                          # Time step, in seconds
    u  = [0, 0, 0]                    # No control input

    x₁ = rk4(x₀, J, u, t, dt)         # Generate next value 

    # Ensure quaternion is unit 
    @test norm(x₁[7:10]) ≈ 1.0 
    
    # Ensure that something actually changed 
    @test x₁ ≉ x₀

end

@testset "Generate an uncontrolled orbit" begin
    x₀ = initialize_orbit(; a = 1) 
    J  = [0.3 0 0; 0 0.3 0; 0 0 0.3]  # Arbitrary inertia matrix for the Satellite 
    t  = Epoch(2020, 11, 30)          # Starting time is Nov 30, 2020
    dt = 1.0                          # Time step, in seconds
    u  = [0, 0, 0]                    # No control input

    N = 5700    # Number of time steps. Should be somewhat around one full orbit
    r_hist = zeros(N, 3)
    nr_hist = zeros(N)
    r_hist[1, :] .= x₀[1:3]
    nr_hist[1]    = norm(x₀[1:3])
    x = x₀
    for i = 1:N - 1
        x = rk4(x, J, u, t, dt)
        t += dt                       # Don't forget to update time (not that it really matters...)
        r_hist[i + 1, :] .= x[1:3]
        nr_hist[i + 1]    = norm(x[1:3])
    end

    plot(r_hist, label = ["x" "y" "z"], title = "Satellite Orbit Test", xlabel = "Time", ylabel = "Satellite Position (m)")
    display(plot!(nr_hist, ls = :dash, label = "Magnitude"))

    @test abs(nr_hist[1] - nr_hist[end]) < 0.001 * nr_hist[1]   # Make sure it ends up in about the same spot
end

@testset "Set a satellite spinning" begin
    x₀ = initialize_orbit() 
    x₀[11:13] .= zeros(3)  # No angular velocity 

    J  = [0.3 0 0; 0 0.3 0; 0 0 0.3]  # Arbitrary inertia matrix for the Satellite 
    t  = Epoch(2020, 11, 30)          # Starting time is Nov 30, 2020
    dt = 0.5                          # Time step, in seconds
    u  = [0.01, 0, 0]                    # Set it spinning

    N = 200
    q_hist = zeros(N, 4)
    q_hist[1, :] .= x₀[7:10]
    x = x₀
    for i = 1:N - 1
        x = rk4(x, J, u, t, dt)
        t += dt                      # Don't forget to update time (not that it really matters...)
        q_hist[i + 1, :] .= x[7:10]
    end

    display(plot(q_hist, title = "Satellite Attitude Test", xlabel = "Time", ylabel = "Quaternion magic"))
end

@testset "Energy conservation" begin 
    function energy(x, J; μ = 3.9860044188e14)
        r, v, q, ω = x[1:3], x[4:6], x[7:10], x[11:13]

        Eₜ = (0.5) * (v' * v)  - (μ / norm(r))   # Specific orbital energy, Vis-Viva equation 
        Eᵣ = 0.5 * ω' * J * ω
        return Eₜ + Eᵣ
    end

    x₀ = initialize_orbit(; a = 1) 
    J  = [0.3 0 0; 0 0.3 0; 0 0 0.3]  # Arbitrary inertia matrix for the Satellite 
    t  = Epoch(2020, 6, 3)            # Starting time 
    dt = 1.0                          # Time step, in seconds
    u  = [0, 0, 0]                    # No control input

    N  = 11000   # Something like 2ish orbits 
    e_hist = zeros(N)
    e_hist[1] = energy(x₀, J)
    x = x₀ 
    for i = 1:N - 1
        x = rk4(x, J, u, t, dt)
        t += dt                       # Don't forget to update time (not that it really matters...)
        e_hist[i + 1] = energy(x, J)
    end

    display(plot(e_hist, title ="Energy"))
end

@testset "MEKF/DeTumbling" begin 
    function control_law(ω, b)
        b̂ = b / norm(b)
        k = 7e-4
        M = -k * (I(3) - b̂*b̂')*ω
        m = 1 / (dot(b, b)) * cross(b, M)
        return cross(m, b) 
    end

    data = my_sim_kf(control_law, 5)
    display(plot(
        data, 
        title="MEKF/DeTumbling", 
        xlabel="Time (s)", 
        ylabel="Angular Velocity (rad/s)", 
        labels=["ω" "s" "v1" "v2" "v3" "s'" "v1'" "v2'" "v3'"],
        linecolor=[:red :blue :green :purple :orange :blue :green :purple :orange]
    ))

end


@testset "DeTumbling" begin 
    function control_law(ω, b)
        b̂ = b / norm(b)
        k = 7e-4
        M = -k * (I(3) - b̂*b̂')*ω
        m = 1 / (dot(b, b)) * cross(b, M)
        return cross(m, b) 
    end

    data = my_sim(control_law)
    display(plot(data, title="DeTumbling", xlabel="Time (s)", ylabel="Angular Velocity (rad/s)", labels=["ω1" "ω2" "ω3" "ω"]))

end


@testset "DeTumbleIO" begin 
    println("DeTumbleIO")
    inp = Pipe()
    out = Pipe()
    err = Pipe()
    proc = run(Cmd(`python3 -u main.py`, dir = "/home/thetazero/Documents/pycubed/software_example_beepsat/state_machine/build/" ), inp, out, err, wait = false)
    close(out.in)
    close(err.in)
    Base.start_reading(out.out)
    Base.start_reading(err.out)

    println(readline(proc.out))
    
    i=0

    function control_law(ω, b)
        println("wrote to sensors")
        control = [0,0,0]
        while true
            write(inp, ">>>ω"*string(ω)*"\n")
            write(inp, ">>>b"*string(b)*"\n")
            write(inp, ">>>?\n")
            input = readline(proc.out)
            println(input)
            if length(input) >= 4 && input[1:3] == ">>>"
                if input[4] == 'M'
                    control = (input[5:length(input)])
                    control = JSON.parse(control)
                end
                break
            end
        end
        i+=1
        println(i)
        return control
    end

    data = my_sim(control_law)
    display(plot(data, title="DeTumbleIO", xlabel="Time (s)", ylabel="Angular Velocity (rad/s)", labels=["ω1" "ω2" "ω3" "ω"]))

end
