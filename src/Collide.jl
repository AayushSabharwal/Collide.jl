module Collide

using ModelingToolkit, DifferentialEquations
using StaticArrays
using LinearAlgebra
using Reexport

@reexport import DifferentialEquations: step!
@reexport using PrimitiveCollisions

@variables t
const D = Differential(t)

export Entity, World, get_system
include("entity.jl")
export Simulation
include("simulate.jl")

using PrecompileTools: PrecompileTools

PrecompileTools.@compile_workload begin
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
    haskey(w, :a)
    w[:a]
    delete!(w, :a)
    push!(w, e)

    sim = Simulation(w, [0.0, 0.0])
    step!(sim, 5.0, true)
end
end
