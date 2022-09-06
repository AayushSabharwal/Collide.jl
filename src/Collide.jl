module Collide

using ModelingToolkit, DifferentialEquations
using PrimitiveCollisions
using StaticArrays
using AbstractTrees
using LinearAlgebra

@variables t
const D = Differential(t)

include("constraints.jl")
include("entity.jl")
include("simulate.jl")
end
