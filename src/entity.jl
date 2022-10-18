"""
    function Entity(; kwargs...)

Create an `Entity` representing a body in the simulation.

# Keywords
- `name::Symbol`: The name for the entity. Necessary, and must be distinct from all other
  entities in the [`World`](@ref).
- `shape::AbstractShape{2} = nothing`: The shape of the entity, specified as a 2D
  primitive from `PrimitiveCollisions.jl`. Can also be `nothing`, which removes the entity
  from the simulation.
- `position::SVector{2,Float64} = @SVector[0.0, 0.0]`: The initial position of the entity.
- `velocity::SVector{2,Float64} = @SVector[0.0, 0.0]`: The initial velocity of the entity.
- `rotation::Float64 = 0.0`: The initial rotation of the entity. Specified in radians in
  an anticlockwise direction from the positive x-axis.
- `angular_velocity::Float64 = 0.0`: The initial angular velocity of the entity. Specified
  in radians per second, in the same direction as `rotation`.
- `mass::Float64 = 1.0`: The mass of the entity. Can be `Inf` to make the body immovable.
- `inertia::Float64 = 1.0`: The moment of inertia ("rotational mass") of the entity. Can
  be `Inf` to make the body impossible to rotate.
- `linear_drag::Float64 = 0.1`: The linear drag coefficient. Linear drag force is calculated
  as the product of this coefficient and the negative of the velocity of the entity.
- `angular_drag::Float64 = 0.1`: The angular drag coefficient. Angular drag torque is
  calculated as the product of this coefficient and the negative of the angular velocity
  of the entity.
- `bounce::Float64 = 1.0`: The coefficient of restitution used during collision resolution
  is taken as the minimum of the `bounce` values of the two involved bodies. Typically in
  the range `[0.0, 1.0]`. 1.0 is a perfectly elastic collision, 0.0 is a perfectly
  inelastic collision.

# Symbolic variables
In the `ODESystem` returned by [`get_system`](@ref) the variables/parameters for an entity
are referred to by sightly different names, prefixed by the name of the entity followed by
an underscore (`_`). Following are the names (`\$name_` is the prefix):
- `\$name_pos(t)[1:2]`: Position
- `\$name_vel(t)[1:2]`: velocity
- `\$name_acc(t)[1:2]`: Acceleration
- `\$name_θ(t)`: Rotation
- `\$name_ω(t)`: Angular velocity
- `\$name_α(t)`: Angular acceleration
- `\$name_m`: Mass
- `\$name_I`: Moment of inertia
- `\$name_μ`: Linear drag coefficient
- `\$name_η`: Angular drag coefficient
"""
Base.@kwdef struct Entity{S<:Union{AbstractShape{2},Nothing}}
    name::Symbol
    shape::S = nothing
    position::SVector{2,Float64} = zero(SVector{2,Float64})
    velocity::SVector{2,Float64} = zero(SVector{2,Float64})
    rotation::Float64 = 0.0
    angular_velocity::Float64 = 0.0
    # TODO automatic mass and inertia calculation
    mass::Float64 = 1.0
    inertia::Float64 = 1.0
    linear_drag::Float64 = 0.1
    angular_drag::Float64 = 0.1
    bounce::Float64 = 1.0
    function Entity(
        name::Symbol,
        shape::S,
        position,
        velocity,
        rotation,
        angular_velocity,
        mass,
        inertia,
        linear_drag,
        angular_drag,
        bounce,
    ) where {S}
        if name == :world
            error("Name cannot be :world")
        end

        return new{S}(
            name,
            shape,
            position,
            velocity,
            rotation,
            angular_velocity,
            mass,
            inertia,
            linear_drag,
            angular_drag,
            bounce,
        )
    end
end

"""
    function Entity(template::Entity; kwargs...)

Create an copy of an existing entity `template`, with some values modified. Keyword
arguments are the same as in the keyword constructor, and `name` must be specified
regardless.
"""
function Entity(
    template::Entity;
    name::Symbol,
    shape = template.shape,
    position = template.position,
    velocity = template.velocity,
    rotation = template.rotation,
    angular_velocity = template.angular_velocity,
    mass = template.mass,
    inertia = template.inertia,
    linear_drag = template.linear_drag,
    angular_drag = template.angular_drag,
    bounce = template.bounce,
)
    return Entity(;
        name,
        shape,
        position,
        velocity,
        rotation,
        angular_velocity,
        mass,
        inertia,
        linear_drag,
        angular_drag,
        bounce,
    )
end

"""
    function get_symbols(ent::Entity)

Get the symbolic variables and parameters for the given entity, in a `NamedTuple`.
"""
function get_symbols(ent::Entity)
    arrsymbols = [:pos, :vel, :acc]
    arrvalues = [ent.position, ent.velocity, [0.0, 0.0]]
    sclsymbols = [:θ, :ω, :α]
    sclvalues = [ent.rotation, ent.angular_velocity, 0.0]
    psymbols = [:m, :I, :μ, :η]
    pvalues = [ent.mass, ent.inertia, ent.linear_drag, ent.angular_drag]

    sts = []
    for (sym, val) in zip(arrsymbols, arrvalues)
        nm = Symbol(ent.name, :_, sym)
        push!(sts, (@variables ($nm)(t)[1:2] = val)[1])
    end

    for (sym, val) in zip(sclsymbols, sclvalues)
        nm = Symbol(ent.name, :_, sym)
        push!(sts, (@variables ($nm)(t) = val)[1])
    end

    prs = []
    for (sym, val) in zip(psymbols, pvalues)
        nm = Symbol(ent.name, :_, sym)
        push!(prs, (@parameters $nm = val)[1])
    end

    sts = NamedTuple{Tuple(vcat(arrsymbols, sclsymbols))}(Tuple(sts))
    prs = NamedTuple{Tuple(psymbols)}(Tuple(prs))
    return sts, prs
end

"""
    function get_self_equations(ent::Entity{S}, sts, prs, gravity)

Given an entity, its symbolic states and parameters, and the force of gravity,
returns the equations of motion for the entity.
"""
function get_self_equations(::Entity{S}, sts, prs, gravity) where {S<:AbstractShape{2}}
    return [
        collect(D.(sts.pos) .~ sts.vel)
        collect(D.(sts.vel) .~ sts.acc)
        collect(sts.acc .~ gravity .- prs.μ .* sts.vel ./ prs.m)
        D(sts.θ) ~ sts.ω
        D(sts.ω) ~ sts.α
        sts.α ~ (-prs.η * sts.ω / prs.I)
    ]
end

abstract type AbstractConstraint end

"""
    function World(name::Symbol)

Creates a world containing all [`Entity`](@ref)s to be simulated. The system returned from
[`get_system`](@ref) is named `name`. Entities can be added to a world using
[`Base.push!`](@ref). [`Base.getindex`](@ref) returns an [`Entity`](@ref) given its name.
[`Base.haskey`](@ref) checks whether an entity with the given name exists in this world.
[`Base.delete!`](@ref) deletes the entity with the given name from the world.
"""
struct World
    name::Symbol
    entities::Dict{Symbol,Entity}
    constraints::Dict{Set{Symbol},AbstractConstraint}
end

World(name::Symbol) = World(name, Dict(), Dict())

"""
    function Base.haskey(w::World, sym::Symbol)

Check whether the [`World`](@ref) contains an entity with name `sym`.
"""
Base.haskey(w::World, sym::Symbol) = haskey(w.entities, sym)

Base.haskey(w::World, a::Symbol, b::Symbol) = haskey(w.constraints, Set((a, b)))

"""
    function Base.getindex(w::World, sym::Symbol)

Return the [`Entity`](@ref) with name `sym` in [`World`](@ref).
"""
Base.getindex(w::World, sym::Symbol) = w.entities[sym]

Base.getindex(w::World, a::Symbol, b::Symbol) = w.constraints[Set((a, b))]

"""
    function Base.push!(w::World, e::Entity)

Add the given [`Entity`](@ref) to the given [`World`](@ref).
"""
function Base.push!(w::World, e::Entity)
    haskey(w, e.name) && error("Entity with name $(e.name) already exists")
    w.entities[e.name] = e
    return w
end

function Base.setindex!(w::World, c::AbstractConstraint, a::Symbol, b::Symbol)
    key = Set((a, b))
    haskey(w.constraints, key) || (w.constraints[key] = AbstractConstraint[])
    push!(w.constraints[key], c)
    return c
end

"""
    function Base.delete!(w::World, sym::Symbol)

Remove the [`Entity`](@ref) with the name `sym` from the given [`World`](@ref).
"""
function Base.delete!(w::World, sym::Symbol)
    delete!(w.entities, sym)
    delete!.((w.constraints,), findall(k -> sym in k, w.constraints))
    return w
end

function Base.delete!(w::World, a::Symbol, b::Symbol)
    delete!(w.constraints, Set((a, b)))
    return w
end

"""
    function get_system(w::World, gravity = [0.0, -9.81])

Return the `ODESystem` for the given [`World`]. Uses the provided value for the
acceleration due to gravity, exposed in the `ODESystem` as `g`.
"""
function get_system(w::World, gravity = [0.0, -9.81])
    @parameters g[1:2] = gravity

    sts = Dict{Symbol,NamedTuple}()
    prs = Dict{Symbol,NamedTuple}()
    eqs = Equation[]
    for ename in keys(w.entities)
        s, p = get_symbols(w[ename])
        sts[ename] = s
        prs[ename] = p
        append!(eqs, get_self_equations(w[ename], s, p, g))
    end

    nested_values(x) = reduce(vcat, collect.(reduce(vcat, collect.(values.(values(x))))))

    sys = ODESystem(
        eqs, t, nested_values(sts), vcat(nested_values(prs), collect(g)); name = w.name
    )
    return sys
end
