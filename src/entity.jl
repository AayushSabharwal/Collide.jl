@kwdef struct Entity{F,S<:Shape{F}}
    name::Symbol
    shape::S
    pos::Point{F} = zero(Point{F})
    vel::Point{F} = zero(Point{F})
    acc::Point{F} = zero(Point{F})
    rot::F = zero(F)
    rvel::F = zero(F)
    racc::F = zero(F)
    m::F = one(F)
    I::F = one(F)
    lindrag::F = F(0.1)
    angdrag::F = F(0.1)
    spring_constant::F = F(100)
    damping_coeff::F = F(10)
end

Entity{F}(; shape, kwargs...) where {F} = Entity{F,typeof(shape)}(; shape, kwargs...)

function get_system(entity::Entity{Fl}, gravity) where {Fl}
    @variables pos(t)[1:2] = entity.pos vel(t)[1:2] = entity.vel acc(t)[1:2] = entity.acc
    @variables θ(t) = entity.rot ω(t) = entity.rvel α(t) = entity.racc
    @variables F(t)[1:2] = [zero(Fl), zero(Fl)] τ(t) = zero(Fl)
    @parameters m = entity.m I = entity.I μ = entity.lindrag η = entity.angdrag
    @parameters k = entity.spring_constant c = entity.damping_coeff


    ODESystem(Equation[], t, [pos..., vel..., acc..., F..., θ, ω, α, τ], [m, I, μ, η, k, c]; name = entity.name)
end

function get_equations(e_sys, gravity)
    return [
        collect(D.(e_sys.pos) .~ e_sys.vel)
        collect(D.(e_sys.vel) .~ e_sys.acc)
        collect(e_sys.acc .~ gravity .- e_sys.μ .* e_sys.vel ./ e_sys.m .+ e_sys.F ./ e_sys.m)
        D(e_sys.θ) ~ e_sys.ω
        D(e_sys.ω) ~ e_sys.α
        e_sys.α ~ -e_sys.η * e_sys.ω / e_sys.I + e_sys.τ / e_sys.I
    ]
end

# i j k
# 0 0 ω
# x y 0
function reverse_transform(point, pos, θ)
    s = sin(θ)
    c = cos(θ)
    rot = SMatrix{2,2}(c, s, -s, c)
    return rot * SVector{2}(point) .+ pos
end

function get_collision_system(ea::Entity{Fl}, eb::Entity{Fl}, ea_sys, eb_sys) where {Fl}

    cur_collision = collision(ea.shape, eb.shape, State(eb.pos .- ea.pos, eb.rot - ea.rot))
    i_pa = reverse_transform(cur_collision[2:3], ea.pos, ea.rot) .- ea.pos
    i_pb = reverse_transform(cur_collision[4:5], ea.pos, ea.rot) .- eb.pos
    cnorm = (i_pb .+ eb.pos .- i_pa .- ea.pos) ./ cur_collision[1]
    i_va = ea.vel .+ SVector{2}(-i_pa[2], i_pa[1]) .* ea.rvel
    i_vb = eb.vel .+ SVector{2}(-i_pb[2], i_pb[1]) .* eb.rvel
    @variables result(t)[1:5] = cur_collision dist(t) = cur_collision[1]
    @variables pa(t)[1:2] = i_pa pb(t)[1:2] = i_pb
    @variables va(t)[1:2] = i_va vb(t)[1:2] = i_vb
    @variables n(t)[1:2] = cnorm
    @variables Fa_raw(t)[1:2] = zero(SVector{2,Fl}) Fb_raw(t)[1:2] = zero(SVector{2,Fl})
    @variables Fa(t)[1:2] = zero(SVector{2,Fl}) Fb(t)[1:2] = zero(SVector{2,Fl})
    @variables τa(t) = zero(Fl) τb(t) = zero(Fl)

    return ODESystem(
        Equation[],
        t,
        [result..., dist, pa..., pb..., va..., vb..., n..., Fa_raw..., Fb_raw..., Fa..., Fb..., τa, τb],
        [];
        name = Symbol(:col_, ea.name, :_, eb.name)
    )
end

# rotates axes by θ
function rot_mat(θ)
    s, c = sincos(θ)
    return SMatrix{2,2}(c, -s, s, c)
end

function connect!(ea::Entity{Fl}, eb::Entity{Fl}, ea_sys, eb_sys, col_sys) where {Fl}
    return [
        collect(col_sys.result .~ collision(ea.shape, eb.shape, State(rot_mat(ea_sys.θ) * SVector{2}(eb_sys.pos .- ea_sys.pos), eb_sys.θ - ea_sys.θ)))
        col_sys.dist ~ col_sys.result[1]
        collect(col_sys.pa .~ rot_mat(-ea_sys.θ) * collect(col_sys.result[2:3]))
        collect(col_sys.pb .~ rot_mat(-ea_sys.θ) * collect(col_sys.result[4:5]) .+ ea_sys.pos .- eb_sys.pos)
        collect(col_sys.va .~ ea_sys.vel .+ SVector{2}(-col_sys.pa[2], col_sys.pa[1]) .* ea_sys.ω)
        collect(col_sys.vb .~ eb_sys.vel .+ SVector{2}(-col_sys.pb[2], col_sys.pb[1]) .* eb_sys.ω)
        collect(col_sys.n .~ (col_sys.pb .+ eb_sys.pos .- col_sys.pa .- ea_sys.pos) ./ (col_sys.dist))
        collect(col_sys.Fb_raw .~ (col_sys.dist <= zero(Fl)) * (ea_sys.k * -col_sys.dist) * col_sys.n)
        collect(col_sys.Fa_raw .~ (col_sys.dist <= zero(Fl)) * (eb_sys.k * -col_sys.dist) * -col_sys.n)
        collect(col_sys.Fa .~ -col_sys.pa * dot(col_sys.Fa_raw, -col_sys.pa) / dot(col_sys.pa, col_sys.pa))
        collect(col_sys.Fb .~ -col_sys.pb * dot(col_sys.Fb_raw, -col_sys.pb) / dot(col_sys.pb, col_sys.pb))
        col_sys.τa ~ col_sys.pa[1] * col_sys.Fa_raw[2] - col_sys.pa[2] * col_sys.Fa_raw[1]
        col_sys.τb ~ col_sys.pb[1] * col_sys.Fb_raw[2] - col_sys.pb[2] * col_sys.Fb_raw[1]
    ]
end

struct World{E<:Entity}
    entities::Vector{E}
end

function get_system(world::World, gravity)
    @parameters g[1:2] = gravity

    systems = Dict{Symbol,ODESystem}()
    force_rhs = Dict{Symbol,SVector{2,Num}}()
    torque_rhs = Dict{Symbol,Num}()
    eqs = Equation[]
    for e in world.entities
        systems[e.name] = get_system(e, g)
        append!(eqs, get_equations(systems[e.name], gravity))
        force_rhs[e.name] = SVector{2}(zero(Num), zero(Num))
        torque_rhs[e.name] = zero(Num)
    end

    collision_sys = Dict{Tuple{Symbol,Symbol},ODESystem}()
    for (ea, eb) in Iterators.product(world.entities, world.entities)
        if ea.name == eb.name || haskey(collision_sys, (eb.name, ea.name))
            continue
        end

        c_sys = get_collision_system(ea, eb, systems[ea.name], systems[eb.name])
        append!(eqs, connect!(ea, eb, systems[ea.name], systems[eb.name], c_sys))
        collision_sys[(ea.name, eb.name)] = c_sys
        force_rhs[ea.name] = force_rhs[ea.name] .+ c_sys.Fa
        force_rhs[eb.name] = force_rhs[eb.name] .+ c_sys.Fb
        torque_rhs[ea.name] = torque_rhs[ea.name] + c_sys.τa
        torque_rhs[eb.name] = torque_rhs[eb.name] + c_sys.τb
    end

    for e in world.entities
        append!(eqs, collect(systems[e.name].F .~ force_rhs[e.name]))
        push!(eqs, systems[e.name].τ ~ torque_rhs[e.name])
    end
    return compose(ODESystem(eqs, t, [], [g...]; name = :world), vcat(collect(values(systems)), collect(values(collision_sys))))
end
