using Collide
using StaticArrays
using Test

@testset "Collide.jl" begin
    e = Entity(; name = :a, shape = Rect(1.0, 1.0), linear_drag = 0.5)
    e2 = Entity(;
        name = :b,
        shape = Rect(1.0, 1.0),
        position = SVector{2}(2.1, -1.2),
        velocity = SVector{2}(-1.0, 0.0),
    )
    w = World(:w)
    push!(w, e)
    push!(w, e2)
    @test haskey(w, :a)
    @test w[:a] == e
    delete!(w, :a)
    @test !haskey(w, :a)
    push!(w, e)

    sim = Simulation(w, [0.0, 0.0])
    @test length(sim.callback.collision_pairs) == 1
    step!(sim, 5.0, true)
    @test sim.integrator.t == 5.0
end
