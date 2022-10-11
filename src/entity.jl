# TODO nothing entity can only be root, otherwise merge children

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

struct World
    name::Symbol
    entities::Dict{Symbol,Entity}
    constraints::Dict{Set{Symbol},AbstractConstraint}
end

World(name) = World(name, Dict(), Dict())

Base.haskey(w::World, sym::Symbol) = haskey(w.entities, sym)

Base.haskey(w::World, a::Symbol, b::Symbol) = haskey(w.constraints, Set((a, b)))

Base.getindex(w::World, sym::Symbol) = w.entities[sym]

Base.getindex(w::World, a::Symbol, b::Symbol) = w.constraints[Set((a, b))]

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

function Base.delete!(w::World, sym::Symbol)
    delete!(w.entities, sym)
    delete!.((w.constraints,), findall(k -> sym in k, w.constraints))
    return w
end

function Base.delete!(w::World, a::Symbol, b::Symbol)
    delete!(w.constraints, Set((a, b)))
    return w
end

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

    sys = ODESystem(eqs, t, nested_values(sts), vcat(nested_values(prs), collect(g)); name = w.name)
    return sys
end
