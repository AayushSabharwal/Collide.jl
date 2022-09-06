using SciMLBase

struct CallbackData
    world::World
    system::ODESystem
    collision_pairs::Vector{NTuple{2,Symbol}}
    collision_axes::Vector{SVector{2,Float64}}
end

struct Simulation{P<:DEProblem,I<:SciMLBase.DEIntegrator}
    problem::P
    integrator::I
    callback::CallbackData
end

function Simulation(world::World, gravity = [0.0, -9.81])
    sys = structural_simplify(get_system(world, gravity))

    collision_pairs = Set{NTuple{2,Symbol}}()
    # eqs = Equation[]
    # events = Pair{Vector{Equation}}[]
    # states = Num[]
    for n1 in PreOrderDFS(IndexNode(world))
        for n2 in PreOrderDFS(IndexNode(world))
            if n1.index == n2.index ||
                !haskey(world.entities, n1.index) ||
                !haskey(world.entities, n2.index) ||
                isnothing(world[n1.index].shape) ||
                isnothing(world[n2.index].shape) ||
                (n2.index, n1.index) in collision_pairs ||
                (n1.index, n2.index) in collision_pairs
                continue
            end

            p1 = AbstractTrees.parent(n1)
            p2 = AbstractTrees.parent(n2)

            if isnothing(p1)
                p2.index == world.name && continue
                push!(collision_pairs, (n1.index, n2.index))
            end

            if isnothing(p2)
                p1.index == world.name && continue
                push!(collision_pairs, (n1.index, n2.index))
            end

            if p1.index == p2.index && !isnothing(world[world.name].shape) ||
                p1.index == n2.index ||
                p2.index == n1.index
                continue
            end

            push!(collision_pairs, (n1.index, n2.index))
            # dist_name = Symbol("sep_", n1.index, "_", n2.index)
            # ax_name = Symbol("dir_", n1.index, "_", n2.index)
            # pa_name = Symbol("closest_", n1.index, "_", n2.index)
            # pb_name = Symbol("closest_", n2.index, "_", n1.index)
            # temp_state = State(
            #     world[n2.index].position - world[n1.index].position,
            #     world[n2.index].rotation - world[n1.index].rotation,
            # )
            # coldata = check_collision(
            #     world[n1.index].shape, world[n2.index].shape, temp_state
            # )
            # pa, pb = closest_pair(world[n1.index].shape, world[n2.index].shape, temp_state)
            # vars = @variables $dist_name(t) = coldata.separation $ax_name(t)[1:2] =
            #     coldata.direction $pa_name(t)[1:2] = pa $pb_name(t)[1:2] = pb

            # sym_state = State(
            #     SVector{2}(
            #         collect(
            #             get_system_property(sys, n2.index, :pos) .-
            #             get_system_property(sys, n1.index, :pos),
            #         ),
            #     ),
            #     get_system_property(sys, n2.index, :θ) -
            #     get_system_property(sys, n1.index, :θ),
            # )
            # sym_col = check_collision(
            #     world[n1.index].shape, world[n2.index].shape, sym_state
            # )
            # sym_pa, sym_pb = closest_pair(
            #     world[n1.index].shape, world[n2.index].shape, sym_state
            # )
            # append!(
            #     eqs,
            #     [
            #         vars[1] ~ sym_col.separation
            #         collect(vars[2] .~ sym_col.direction)
            #         collect(vars[3] .~ sym_pa)
            #         collect(vars[4] .~ sym_pb)
            #     ],
            # )
            # append!(states, [vars[1], vars[2]..., vars[3]..., vars[4]...])
            # push!(
            #     events,
            #     [vars[1] ~ 0.0] => (
            #         affect!,
            #         [
            #             get_system_property(sys, n1.index, :pos) => :posa
            #             get_system_property(sys, n2.index, :pos) => :posb
            #             get_system_property(sys, n1.index, :vel) => :ua
            #             get_system_property(sys, n2.index, :vel) => :ub
            #             get_system_property(sys, n1.index, :θ) => :rota
            #             get_system_property(sys, n2.index, :θ) => :rotb
            #             vars[2] => :cnorm
            #             vars[3] => :cpa
            #             vars[4] => :cpb
            #         ],
            #         [
            #             get_system_property(sys, n1.index, :m) => :ma
            #             get_system_property(sys, n2.index, :m) => :mb
            #         ],
            #         (n1.index, n2.index),
            #     ),
            # )
        end
    end
    # collision_sys = ODESystem(eqs, t, states, []; name = :collision, continuous_events = events)
    # sys = structural_simplify(compose(sys, collision_sys))
    problem = ODAEProblem(sys, [], (0.0, Inf))
    cbdata = CallbackData(
        world,
        sys,
        collect(collision_pairs),
        zeros(SVector{2,Float64}, length(collision_pairs)),
    )
    if length(collision_pairs) == 0
        integrator = init(problem, Tsit5(); save_everystep = false)
    else
        cb = VectorContinuousCallback(
            cbdata, cbdata, length(collision_pairs); rootfind = SciMLBase.LeftRootFind
        )
        integrator = init(problem, Tsit5(); save_everystep = false, callback = cb)
    end
    return Simulation(problem, integrator, cbdata)
end

get_item_symbol(sysname, itemname) = Symbol("$(sysname)_$(itemname)(t)")

function get_array_item_symbol(sysname, itemname, i)
    return Symbol("getindex(", get_item_symbol(sysname, itemname), ", $i)")
end
function get_array_item_indices(sysname, itemname, integrator)
    return SciMLBase.sym_to_index.(
        get_array_item_symbol.(sysname, itemname, (1, 2)), (integrator,)
    )
end

function (cbdata::CallbackData)(out, u, t, integrator)
    for i in 1:length(out)
        body_a, body_b = cbdata.collision_pairs[i]
        shape_a = cbdata.world[body_a].shape
        shape_b = cbdata.world[body_b].shape
        i_pos_a = get_array_item_indices(body_a, "pos", integrator)
        i_pos_b = get_array_item_indices(body_b, "pos", integrator)
        pos_a = getindex.((u,), i_pos_a)
        pos_b = getindex.((u,), i_pos_b)
        i_rot_a = SciMLBase.sym_to_index(get_item_symbol(body_a, "θ"), integrator)
        i_rot_b = SciMLBase.sym_to_index(get_item_symbol(body_b, "θ"), integrator)
        rot_a = getindex(u, i_rot_a)
        rot_b = getindex(u, i_rot_b)
        state = State(SVector{2}(pos_b .- pos_a), rot_b - rot_a)
        coldata = check_collision(shape_a, shape_b, state)
        out[i] = coldata.separation
        cbdata.collision_axes[i] = coldata.direction
        # @show t pos_a pos_b state coldata
    end
end

get_system_property(system, obj, property) = getproperty(system, Symbol(obj, "_", property))

function (cbdata::CallbackData)(integrator, eidx)
    ba, bb = cbdata.collision_pairs[eidx]
    ea = cbdata.world[ba]
    eb = cbdata.world[bb]
    symua = get_system_property(cbdata.system, ba, "vel")
    symub = get_system_property(cbdata.system, bb, "vel")
    println(symua, symub)
    ua = integrator[symua]
    ub = integrator[symub]
    ma = ea.mass
    mb = eb.mass
    Ia = ea.inertia
    Ib = eb.inertia
    posa = integrator[get_system_property(cbdata.system, ba, "pos")]
    posb = integrator[get_system_property(cbdata.system, bb, "pos")]
    rota = integrator[get_system_property(cbdata.system, ba, "θ")]
    rotb = integrator[get_system_property(cbdata.system, bb, "θ")]
    ang_ua = integrator[get_system_property(cbdata.system, ba, "ω")]
    ang_ub = integrator[get_system_property(cbdata.system, bb, "ω")]
    e = min(ea.bounce, eb.bounce)
    cnorm = cbdata.collision_axes[eidx]
    s, c = sincos(rota)
    rot_mat = SMatrix{2,2}(c, s, -s, c)
    cnorm = rot_mat * cnorm
    state = State(SVector{2}(posb - posa), rotb - rota)
    pa, pb = closest_pair(ea.shape, eb.shape, state)
    @show pa, pb
    pa = rot_mat * pa + posa
    pb = rot_mat * pb + posa
    @show norm(pa - pb)
    col = check_collision(ea.shape, eb.shape, state)
    @show col
    r1 = SVector{3}((pa - posa)..., 0.0)
    r2 = SVector{3}((pb - posb)..., 0.0)

    vc = (
        ua + cross(SVector{3}(0.0, 0.0, ang_ua), r1)[1:2] - ub -
        cross(SVector{3}(0.0, 0.0, ang_ub), r2)[1:2]
    )

    J =
        -(1.0 + e) * dot(vc, cnorm) / (
            1.0 / ma +
            1.0 / mb +
            dot(
                cross(cross(r1, SVector{3}(cnorm..., 0.0)), r1) / Ia +
                cross(cross(r2, SVector{3}(cnorm..., 0.0)), r2) / Ib,
                SVector{3}(cnorm..., 0.0),
            )
        )
    va = ua + J * cnorm / ma
    vb = ub - J * cnorm / mb

    vrota = ang_ua + cross(r1, SVector{3}(cnorm..., 0.0))[3] * J / Ia
    vrotb = ang_ub - cross(r2, SVector{3}(cnorm..., 0.0))[3] * J / Ib

    # TODO use mass
    symva = get_array_item_symbol.(ba, "vel", (1, 2))
    symvb = get_array_item_symbol.(bb, "vel", (1, 2))

    set_u!.((integrator,), symva, va)
    set_u!.((integrator,), symvb, vb)

    set_u!.((integrator,), get_item_symbol(ba, "ω"), vrota)
    set_u!.((integrator,), get_item_symbol(bb, "ω"), vrotb)

    return @show integrator.t ba bb cnorm ua ub va vb pa pb posa posb rota rotb vrota vrotb J
end

# function affect!(integrator, u, p, ctx)
#     ba, bb = ctx
#     ua, ub = u.ua, u.ub
#     posa, posb = u.posa, u.posb
#     rota, rotb = u.rota, u.rotb
#     cnorm = u.cnorm
#     ctang = SMatrix{2,2}(0.0, -1.0, 1.0, 0.0) * cnorm
# 
#     state = State(SVector{2}(posb - posa), rotb - rota)
#     pa, pb = u.cpa, u.cpb
#     atang = dot(ua, ctang)
#     btang = dot(ub, ctang)
#     anorm = dot(ua, cnorm)
#     bnorm = dot(ub, cnorm)
# 
#     # TODO use masses for conservation of momentum
#     # TODO restitution
#     va = cnorm * bnorm + ctang * atang
#     vb = cnorm * anorm + ctang * btang
# 
#     # TODO use mass
#     ΔLa = cross(SVector{3}((pa - posa)..., 0.0), SVector{3}((va - ua)..., 0.0))[3]
#     ΔLb = cross(SVector{3}((pb - posb)..., 0.0), SVector{3}((vb - ub)..., 0.0))[3]
#     symva = get_array_item_symbol.(ba, "vel", (1, 2))
#     symvb = get_array_item_symbol.(bb, "vel", (1, 2))
# 
#     set_u!.((integrator,), symva, va)
#     set_u!.((integrator,), symvb, vb)
# 
#     set_u!.((integrator,), get_item_symbol(ba, "ω"), u.ωa + ΔLa)
#     set_u!.((integrator,), get_item_symbol(bb, "ω"), u.ωb + ΔLb)
# 
#     return println(
#         ba,
#         bb,
#         " ",
#         cnorm,
#         ctang,
#         anorm,
#         " ",
#         atang,
#         " ",
#         bnorm,
#         " ",
#         btang,
#         ua,
#         ub,
#         va,
#         vb,
#         ΔLa,
#         " ",
#         ΔLb,
#     )
# end

step!(sim::Simulation, args...) = DifferentialEquations.step!(sim.integrator, args...)
