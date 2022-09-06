abstract type AbstractConstraint end

struct NoConstraint <: AbstractConstraint end

struct StaticConstraint <: AbstractConstraint end

struct RodConstraint <: AbstractConstraint
    parent_end::SVector{2,Float64}
    self_end::SVector{2,Float64}
end
