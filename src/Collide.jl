module Collide

using ModelingToolkit, DifferentialEquations
using StaticArrays
using LinearAlgebra
# using Reexport

# @reexport import DifferentialEquations: step!

@variables t
const D = Differential(t)

export Point
const Point{F} = SVector{2,F}

export Shape, Circle, Capsule, collision, State
include("collision.jl")
export Entity, World, get_system, get_collision_system
include("entity.jl")

include("plotting.jl")
# export Simulation
# include("simulate.jl")

# using SnoopPrecompile: SnoopPrecompile

end
