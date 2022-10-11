using SciMLBase

const IndexesTupleType = NamedTuple{(:pos, :vel, :θ, :ω, :m, :I, :μ, :η, :entity), Tuple{SVector{2, Int}, SVector{2, Int}, Int, Int, Int, Int, Int, Int, Int}};

struct CallbackData{S}
    system::ODESystem
    collision_pairs::Vector{NTuple{2,Symbol}}
    collision_axes::Vector{SVector{2,Float64}}
    symbols::Dict{Symbol,IndexesTupleType}
    shapes::Vector{S}
    bounces::Vector{Float64}
end

struct Simulation{P<:DEProblem,I<:SciMLBase.DEIntegrator,C<:CallbackData}
    problem::P
    integrator::I
    callback::C
end

function Simulation(world::World, gravity = [0.0, -9.81])
    sys = structural_simplify(get_system(world, gravity))

    collision_pairs = Set{NTuple{2,Symbol}}()
    shape_type = Union{}
    # collision_checks = Function[]
    for a in keys(world.entities)
        for b in keys(world.entities)
            if a == b ||
               isnothing(world[a].shape) ||
               isnothing(world[b].shape) ||
               (a, b) in collision_pairs ||
               (b, a) in collision_pairs
                continue
            end

            push!(collision_pairs, (a, b))
        end

        shape_type = Union{shape_type,typeof(world[a].shape)}
    end

    problem = ODEProblem(sys, [], (0.0, Inf))
    cbdata = CallbackData{shape_type}(
        sys,
        collect(collision_pairs),
        zeros(SVector{2,Float64}, length(collision_pairs)),
        Dict{Symbol,IndexesTupleType}(),
        shape_type[],
        Float64[],
    )
    if length(collision_pairs) == 0
        integrator = init(problem, AutoTsit5(Rosenbrock23()); save_everystep = false)
    else
        cb = VectorContinuousCallback(
            cbdata,
            cbdata,
            length(collision_pairs);
            rootfind = SciMLBase.LeftRootFind,
        )
        integrator =
            init(problem, AutoTsit5(Rosenbrock23()); save_everystep = false, callback = cb)
    end

    for (i, a) in enumerate(keys(world.entities))
        cbdata.symbols[a] = (
            pos = SVector{2}(SciMLBase.sym_to_index.(Symbol.("getindex($(a)_pos(t), " .* ("1", "2") .* ")"), (integrator, ))),
            vel = SVector{2}(SciMLBase.sym_to_index.(Symbol.("getindex($(a)_vel(t), " .* ("1", "2") .* ")"), (integrator, ))),
            θ = SciMLBase.sym_to_index(Symbol("$(a)_θ(t)"), integrator),
            ω = SciMLBase.sym_to_index(Symbol("$(a)_ω(t)"), integrator),
            m = findfirst(==(Symbol("$(a)_m")), SciMLBase.getparamsyms(integrator)),
            I = findfirst(==(Symbol("$(a)_I")), SciMLBase.getparamsyms(integrator)),
            μ = findfirst(==(Symbol("$(a)_μ")), SciMLBase.getparamsyms(integrator)),
            η = findfirst(==(Symbol("$(a)_η")), SciMLBase.getparamsyms(integrator)),
            entity = i
        )
        push!(cbdata.shapes, world[a].shape)
        push!(cbdata.bounces, world[a].bounce)
    end

    return Simulation(problem, integrator, cbdata)
end

function (cbdata::CallbackData)(out, u, t, integrator)
    for i in eachindex(out)
        body_a, body_b = cbdata.collision_pairs[i]
        idxs_a = cbdata.symbols[body_a]
        idxs_b = cbdata.symbols[body_b]
        pos_a = getindex.((u,), idxs_a.pos)
        pos_b = getindex.((u,), idxs_b.pos)
        rot_a = getindex(u, idxs_a.θ)
        rot_b = getindex(u, idxs_b.θ)
        s, c = sincos(rot_a)
        rot_mat = SMatrix{2,2}(c, -s, s, c)
        state = State(rot_mat * SVector{2}(pos_b .- pos_a), rot_b - rot_a)
        coldata::CollisionData{Float64} = check_collision(cbdata.shapes[idxs_a.entity], cbdata.shapes[idxs_b.entity], state)
        out[i] = coldata.separation
        cbdata.collision_axes[i] = SMatrix{2,2}(c, s, -s, c) * coldata.direction
    end
end

function (cbdata::CallbackData)(integrator, eidx)
    ba, bb = cbdata.collision_pairs[eidx]
    idxs_a = cbdata.symbols[ba]
    idxs_b = cbdata.symbols[bb]
    
    ua = getindex.((integrator,), idxs_a.vel)
    ub = getindex.((integrator,), idxs_b.vel)
    ma = integrator.p[idxs_a.m]
    mb = integrator.p[idxs_b.m]
    Ia = integrator.p[idxs_a.I]
    Ib = integrator.p[idxs_b.I]
    posa = getindex.((integrator,), idxs_a.pos)
    posb = getindex.((integrator,), idxs_b.pos)
    rota = integrator[idxs_a.θ]
    rotb = integrator[idxs_b.θ]
    ang_ua = integrator[idxs_a.ω]
    ang_ub = integrator[idxs_b.ω]
    e = min(cbdata.bounces[idxs_a.entity], cbdata.bounces[idxs_b.entity])
    cnorm = cbdata.collision_axes[eidx]

    s, c = sincos(rota)
    rot_mat = SMatrix{2,2}(c, -s, s, c)
    inv_rot_mat = SMatrix{2,2}(c, s, -s, c)
    state = State(rot_mat * SVector{2}(posb - posa), rotb - rota)
    ps::NTuple{2,SVector{2,Float64}} = closest_pair(cbdata.shapes[idxs_a.entity], cbdata.shapes[idxs_b.entity], state)
    pa = inv_rot_mat * ps[1] + posa
    pb = inv_rot_mat * ps[2] + posa
    r1 = pa - posa
    r2 = pb - posb

    my_cross(k, v::SVector{2}) = SVector{2}(-v[2] * k, v[1] * k)
    my_cross(a::SVector{2}, b::SVector{2}) = a[1] * b[2] - a[2] * b[1]

    vc = (ua + my_cross(ang_ua, r1) - ub - my_cross(ang_ub, r2))

    J =
        -(1.0 + e) * dot(vc, cnorm) / (
            1.0 / ma +
            1.0 / mb +
            dot(
                my_cross(my_cross(r1, cnorm), r1) / Ia +
                my_cross(my_cross(r2, cnorm), r2) / Ib,
                cnorm,
            )
        )
    va = ua + J * cnorm / ma
    vb = ub - J * cnorm / mb

    vrota = ang_ua + my_cross(r1, cnorm) * J / Ia
    vrotb = ang_ub - my_cross(r2, cnorm) * J / Ib

    set_u!.((integrator,), getindex.((cbdata.system.states,), idxs_a.vel), va)
    set_u!.((integrator,), getindex.((cbdata.system.states,), idxs_b.vel), vb)

    set_u!(integrator, cbdata.system.states[idxs_a.ω], vrota)
    set_u!(integrator, cbdata.system.states[idxs_b.ω], vrotb)

    return nothing
end

step!(sim::Simulation, args...) = DifferentialEquations.step!(sim.integrator, args...)
