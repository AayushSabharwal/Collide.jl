using GLMakie
import GeometryBasics
using GLMakie.Makie.ColorSchemes
using Random
using JuliaSimCompiler

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

function animate(w::World, sys, duration; framerate = 30, resolution = (1000, 1000), file="animation.mkv", limits=(-10, 10, -10, 10))
    frames = duration * framerate
    timestep = 1 / framerate

    prob = ODEProblem(sys, [], (0.0, Inf))
    int = Observable(init(prob, Tsit5()))

    fig = Figure(; resolution)
    time = Observable(0.0)
    ax = fig[1, 1] = Axis(fig; limits, title=@lift("t=$($time)"))

    cs = ColorSchemes.seaborn_bright
    for (idx, e) in enumerate(w.entities)
        color = getindex(cs, mod1(idx, length(cs)))

        position = @lift(Tuple($int[getproperty(sys, Symbol(e.name, :₊pos))]))
        if e.shape isa Circle
            meshscatter!(ax, position, marker=glCircle(e.shape), color=color, markersize=1, shading=false)
        else
            rot = @lift($int[getproperty(sys, Symbol(e.name, :₊θ))])
            circ_1_pos = @lift($position .+ (cos($rot), sin($rot)) .* e.shape.half_len)
            circ_2_pos = @lift($position .- (cos($rot), sin($rot)) .* e.shape.half_len)

            meshscatter!(@lift([$circ_1_pos, $circ_2_pos]), marker=glCircle(Circle(e.shape.radius)), color=color, markersize=1, shading=false)
            meshscatter!(ax, position, marker=glRect(e.shape), markersize=1, color=color, rotations=rot, shading=false)
        end
        # force = @lift(Tuple($int[getproperty(sys, Symbol(e.name, :₊F))]))
        # lines!(@lift([$position, $position .+ $force]), color=:red)
    end

    # npos = @lift(Tuple($int[sys.col_b_a₊n]))
    # lines!(@lift([(0.0, 0.0), $npos]), color=:green)

    # meshscatter!(@lift([Tuple($int[sys.col_b_a₊pa] .+ $int[sys.b₊pos]), Tuple($int[sys.col_b_a₊pb] .+ $int[sys.a₊pos])]), marker=glCircle(Circle(0.1)), markersize=1)

    record(fig, file, 1:frames; framerate) do _
        time[] = time[] + timestep
        step!(int[], timestep, true)
        int[] = int[]
    end
end

# function manybody(n=10; gravity=[0.0,0.0], k=1000.0, seed=42)
#     w_top = Entity{Float64}(name=:w_top, pos=SVector{2}(0.0, 10.0), m=Inf, I=Inf, shape=Capsule(10.0, 0.5))
#     w_bottom = Entity{Float64}(name=:w_bottom, pos=SVector{2}(0.0, -10.0), m=Inf, I=Inf, shape=Capsule(10.0, 0.5))
#     w_left = Entity{Float64}(name=:w_left, pos=SVector{2}(-10.0, 0.0), rot=π/2, m=Inf, I=Inf, shape=Capsule(10.0, 0.5))
#     w_right = Entity{Float64}(name=:w_right, pos=SVector{2}(10.0, 0.0), rot=π/2, m=Inf, I=Inf, shape=Capsule(10.0, 0.5))

#     entities = Entity{Float64}[w_top, w_bottom, w_left, w_right]
#     rng = MersenneTwister(seed)
#     while n > 0
#         shape = if rand(rng) < 0.5
#             Circle(rand(rng) * 0.8 + 0.2)
#         else
#             hlen = rand(rng) * 0.8 + 0.2
#             rad = rand(rng) * (hlen - 0.2) + 0.2
#             Capsule(hlen, rad)
#         end


#         pos = @SVector[
#             rand(rng) * 17 - 8.5,
#             rand(rng) * 17 - 8.5,
#         ]

#         rot = rand(rng) * 2π

#         vel_ang = rand(rng) * 2π
#         vel = (rand(rng) * 2.0 + 3.0) * @SVector[cos(vel_ang), sin(vel_ang)]

#         rvel = rand(rng) * 3.0

#         m = π * shape.radius ^ 2
#         if shape isa Capsule
#             m += shape.half_len * shape.radius * 4
#         end
#         I = if shape isa Circle
#             1
#         else
#             m * (shape.half_len ^ 2 + shape.radius ^ 2) / 3
#         end

#         possible = true
#         for e in entities
#             if collision(e.shape, shape, State(rot_mat(e.rot) * (pos .- e.pos), rot - e.rot))[1] <= 0
#                 possible = false
#                 break
#             end
#         end
#         possible || continue


#         push!(entities, Entity{Float64}(; name=Symbol(:shape_, n), shape, pos, rot, vel, rvel, m, I))
#         n -= 1
#     end

#     return World(entities)
# end
