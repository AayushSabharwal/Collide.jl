# TODO nothing entity can only be root, otherwise merge children

Base.@kwdef struct Entity{S<:Union{AbstractShape{2},Nothing},C<:AbstractConstraint}
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
    constraint::C = NoConstraint()
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
        constraint::C,
    ) where {S,C}
        if name == :root
            error("Name cannot be :root")
        end

        return new{S,C}(
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
            constraint,
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
    constraint = template.constraint,
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
        constraint,
    )
end

function get_symbols(ent::Entity, parent_pos, gravity)
    arrsymbols = [:pos, :vel, :acc]
    arrvalues = [parent_pos .+ ent.position, ent.velocity, [0.0, 0.0]]
    sclsymbols = [:θ, :ω, :ɑ]
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

function get_self_equations(ent::Entity{S}, sts, prs, gravity) where {S<:AbstractShape{2}}
    return [
        collect(D.(sts.pos) .~ sts.vel)
        collect(D.(sts.vel) .~ sts.acc)
        collect(sts.acc .~ gravity .- prs.μ .* sts.vel ./ prs.m)
        D(sts.θ) ~ sts.ω
        D(sts.ω) ~ sts.ɑ
        sts.ɑ ~ (-prs.η * sts.ω / prs.I)
    ]
end

function get_self_equations(
    ent::Entity{S,StaticConstraint}, sts, prs, gravity
) where {S<:AbstractShape{2}}
    return [
        collect(sts.pos .~ ent.position)
        collect(sts.vel .~ ent.velocity)
        collect(sts.acc .~ [0.0, 0.0])
        sts.θ ~ ent.rotation
        sts.ω ~ ent.angular_velocity
        sts.ɑ ~ 0.0
    ]
end

function get_relative_equations(ent::Entity{S}, args...) where {S<:AbstractShape{2}}
    return Equation[]
end

function get_relative_equations(
    ent::Entity{S,RodConstraint}, selfsts, selfprs, parsts, parprs, parent_pos, gravity
) where {S<:AbstractShape{2}}
    rodlen = norm(ent.constraint.parent_end .- (ent.position + ent.constraint.self_end))
    s, c = sincos(selfsts.θ)
    self_rot = SMatrix{2,2}(c, s, -s, c)
    s, c = sincos(parsts.θ)
    par_rot = SMatrix{2,2}(c, s, -s, c)
    return [
        0.0 ~
            norm(
                (selfsts.pos .+ self_rot * ent.constraint.self_end) .-
                (parsts.pos .+ par_rot * ent.constraint.parent_end),
            ) - rodlen,
    ]
end

struct World
    name::Symbol
    entities::Dict{Symbol,Entity}
    children::Dict{Symbol,Vector{Symbol}}
    parent::Dict{Symbol,Symbol}
end

function World(; name::Symbol)
    return World(name, Dict(name => Entity(; name)), Dict(name => []), Dict())
end

function World(entity::Entity)
    return World(entity.name, Dict(entity.name => entity), Dict(entity.name => []), Dict())
end

World(w::World) = w

Base.length(w::World) = length(w.entities)

Base.getindex(w::World, idx) = w.entities[idx]

function Base.:+(e::Union{Entity,World}, es::Union{Entity,World}...)
    temp_name = Symbol(rand(UInt16))
    w = (World(e), World.(es)...)
    entities = Dict(
        ∪(getproperty.(w, :entities)..., [temp_name => Entity(; name = temp_name)])
    )
    if length(entities) != sum(length.(getproperty.(w, :entities))) + 1
        error("Duplicate entities")
    end

    children = Dict(
        ∪(getproperty.(w, :children)..., [temp_name => collect(getproperty.(w, :name))])
    )
    parent = Dict(
        ∪(getproperty.(w, :parent)..., collect(getproperty.(w, :name) .=> temp_name))
    )

    return World(temp_name, entities, children, parent)
end

function Base.:∘(e::Entity, w::World)
    if haskey(w.entities, w.name) && !isnothing(w[w.name].shape)
        error("Can only combine when world has no root")
    end
    entities = copy(w.entities)
    delete!(entities, w.name)
    entities[e.name] = e

    children = copy(w.children)
    children[e.name] = w.children[w.name]
    delete!(children, w.name)

    parent = copy(w.parent)
    for ch in children[e.name]
        parent[ch] = e.name
    end
    return World(e.name, entities, children, parent)
end

Base.:∘(e::Entity, w::Entity) = e ∘ World(w)

function get_system(w::World, gravity = [0.0, -9.81])
    @parameters g[1:2] = gravity
    sts = Dict{Symbol,NamedTuple}()
    prs = Dict{Symbol,NamedTuple}()
    eqs = Equation[]
    positions = Dict{Symbol,SVector{2,Float64}}(w.name => w[w.name].position)
    for node in PreOrderDFS(IndexNode(w))
        isroot(node) && continue
        if !isempty(children(node))
            positions[node.index] =
                nodevalue(node).position .+ positions[parentindex(node.tree, node.index)]
        end
        if haskey(node.tree.entities, node.index)
            ent = nodevalue(node)
            parind = parentindex(node.tree, node.index)
            s, p = get_symbols(ent, positions[parind], g)
            sts[node.index] = s
            prs[node.index] = p
            append!(eqs, get_self_equations(ent, s, p, g))
            # append!(
            #     eqs,
            #     get_relative_equations(
            #         ent, s, p, sts[parind], prs[parind], positions[parind], g
            #     ),
            # )
        end
    end
    sys = ODESystem(
        eqs,
        t,
        reduce(vcat, reduce.(vcat, collect.(values.(values(sts))))),
        vcat(reduce(vcat, reduce.(vcat, collect.(values.(values(prs))))), collect(g));
        name = w.name,
    )
    return sys
end

# AbstractTree interface
AbstractTrees.childindices(tree::World, index) = get(tree.children, index, ())
AbstractTrees.ParentLinks(::Type{World}) = AbstractTrees.StoredParents()
AbstractTrees.parentindex(tree::World, index) = get(tree.parent, index, nothing)
AbstractTrees.rootindex(tree::World) = tree.name
function AbstractTrees.ischild(
    n1::IndexNode{World,Symbol}, n2::IndexNode{World,Symbol}; equiv = (===)
)
    return equiv(n1.tree.parent[n1.index], n2.index)
end
AbstractTrees.isroot(x::IndexNode{World,Symbol}) = x.tree.name == x.index
