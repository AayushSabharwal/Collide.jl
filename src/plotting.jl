using GLMakie
import GeometryBasics

function glCircle(c::Circle)
    return GeometryBasics.Circle(GeometryBasics.Point(0.0, 0.0), c.radius)
end

function glRect(c::Capsule)
    return GeometryBasics.Rect(-c.half_len, -c.radius, 2c.half_len, 2c.radius)
    # return GeometryBasics.Polygon([
    #     GeometryBasics.Point(c.half_len, c.radius)
    #     GeometryBasics.Point(-c.half_len, c.radius)
    #     GeometryBasics.Point(-c.half_len, -c.radius)
    #     GeometryBasics.Point(c.half_len, -c.radius)
    # ])
end

function animate(w::World, sys::ODESystem, duration; framerate = 30, resolution = (1000, 1000), file="animation.mkv", limits=(-10, 10, -10, 10))
    frames = duration * framerate
    timestep = 1 / framerate

    prob = ODEProblem(sys, [], (0.0, Inf))
    int = Observable(init(prob, Tsit5()))

    fig = Figure(; resolution)
    time = Observable(0.0)
    ax = fig[1, 1] = Axis(fig; limits, title=@lift("t=$($time)"))

    for e in w.entities
        position = @lift(Tuple($int[getproperty(sys, Symbol(e.name, :₊pos))]))
        if e.shape isa Circle
            meshscatter!(ax, position, marker=glCircle(e.shape), markersize=1)
        else
            rot = @lift($int[getproperty(sys, Symbol(e.name, :₊θ))])
            meshscatter!(ax, position, marker=glRect(e.shape), markersize=1, rotations=rot)
            circ_1_pos = @lift($position .+ (cos($rot), sin($rot)) .* e.shape.half_len)
            circ_2_pos = @lift($position .- (cos($rot), sin($rot)) .* e.shape.half_len)

            meshscatter!(@lift([$circ_1_pos, $circ_2_pos]), marker=glCircle(Circle(e.shape.radius)), markersize=1)
        end
        force = @lift(Tuple($int[getproperty(sys, Symbol(e.name, :₊F))]))
        lines!(@lift([$position, $position .+ $force]), color=:red)
    end

    npos = @lift(Tuple($int[sys.col_b_a₊n]))
    lines!(@lift([(0.0, 0.0), $npos]), color=:green)

    meshscatter!(@lift([Tuple($int[sys.col_b_a₊pa] .+ $int[sys.b₊pos]), Tuple($int[sys.col_b_a₊pb] .+ $int[sys.a₊pos])]), marker=glCircle(Circle(0.1)), markersize=1)

    record(fig, file, 1:frames; framerate) do _
        time[] = time[] + timestep
        step!(int[], timestep, true)
        int[] = int[]
    end
end
