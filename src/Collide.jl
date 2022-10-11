module Collide

using ModelingToolkit, DifferentialEquations
using PrimitiveCollisions
using StaticArrays
using LinearAlgebra

@variables t
const D = Differential(t)

include("entity.jl")
include("simulate.jl")
end
