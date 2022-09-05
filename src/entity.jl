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

    function Entity(name::Symbol, shape::S, args...) where {S}
        if name == :root
            error("Name cannot be :root")
        end

        return new{S}(name, shape, args...)
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

function get_system(ent::Entity{S}, parent_pos, gravity) where {S<:AbstractShape{2}}
    @variables pos(t)[1:2] = (parent_pos .+ ent.position) vel(t)[1:2] = ent.velocity acc(t)[1:2] = [
        0.0, 0.0
    ]
    @variables θ(t) = ent.rotation ω(t) = ent.angular_velocity ɑ(t) = 0.0
    @parameters m = ent.mass I = ent.inertia μ = ent.linear_drag η = ent.angular_drag
    # TODO custom force, torque
    eqs = [
        collect(D.(pos) .~ vel)
        collect(D.(vel) .~ acc)
        collect(acc .~ gravity .- μ .* vel ./ m)
        D(θ) ~ ω
        D(ω) ~ ɑ
        ɑ ~ (-η * ω / I)
    ]

    sys = ODESystem(
        eqs, t, [pos..., vel..., acc..., θ, ω, ɑ], [gravity..., m, I, μ, η]; name = ent.name
    )
    return sys
end

get_system(e::Entity{Nothing}, args...) = ODESystem(Equation[], t; name = e.name)

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
    @parameters g[1:2] = SVector{2}(gravity)
    systems = [get_system(w[w.name], zero(SVector{2,Float64}), g)]
    positions = Dict{Symbol,SVector{2,Float64}}(w.name => w[w.name].position)
    for node in PreOrderDFS(IndexNode(w))
        isroot(node) && continue
        if !isempty(children(node)) && !isroot(node)
            positions[node.index] =
                nodevalue(node).position .+ positions[parentindex(node.tree, node.index)]
        end
        if haskey(node.tree.entities, node.index)
            push!(
                systems,
                get_system(
                    nodevalue(node), positions[parentindex(node.tree, node.index)], g
                ),
            )
        end
    end
    sys = compose(systems...; name = w.name)
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
