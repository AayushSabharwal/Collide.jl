using IfElse

abstract type Shape{F} end

struct Circle{F} <: Shape{F}
    radius::F
end

(c::Circle)(p) = sqrt(dot(p, p)) - c.radius

struct Capsule{F} <: Shape{F}
    half_len::F
    radius::F
end

function (c::Capsule)(p)
    # pa = p - (-c.half_len, 0)
    # ba = (2c.half_len, 0)
    # h = clamp(dot(pa, ba) / dot(ba, ba), 0, 1)
    # diff = pa - h * ba
    # from https://iquilezles.org/articles/distfunctions2d/
    F = eltype(p)
    h = min(one(F), max(zero(F), p[1] / 2c.half_len + F(0.5)))
    diff = SVector{2}(p[1] + (1 - 2h) * c.half_len, p[2])
    return sqrt(dot(diff, diff)) - c.radius
end

struct State{F}
    pos::Point{F}
    rot::F
end

function collision(a::Shape, b::Shape, st::State)
  s, c = sincos(st.rot)
  res = collision(b, a, State(SMatrix{2,2}(c, -s, s, c) * -st.pos, -st.rot))
  rot_mat = SMatrix{2,2}(c, s, -s, c)
  SVector{5}(
    res[1],
    (rot_mat * res[4:5] .+ st.pos)...,
    (rot_mat * res[2:3] .+ st.pos)...,
  )
end

function collision(a::Circle, b::Circle, st::State)
    ba = st.pos
    dist = sqrt(dot(ba, ba))
    return SVector{5}(
        dist - a.radius - b.radius,
        (ba * a.radius / dist)...,
        (st.pos - ba * b.radius / dist)...,
    )
end

function collision(a::Capsule, b::Circle, st::State)
    p = st.pos
    F = eltype(p)
    h = min(one(F), max(zero(F), p[1] / 2a.half_len + F(0.5)))

    pt = SVector{2}((2h - 1) * a.half_len, zero(F))
    dir = p - pt
    dist = sqrt(dot(dir, dir))

    return SVector{5}(
        dist - a.radius - b.radius,
        (pt + dir * a.radius / dist)...,
        (p - dir * b.radius / dist)...,
    )
end

clamp01(x) = min(one(x), max(zero(x), x))

function collision(a::Capsule, b::Capsule, st::State)
    F = typeof(st.rot)
    eps = F(1e-6)
    a0 = SVector{2}(-a.half_len, zero(a.half_len))
    a1 = SVector{2}(a.half_len / 2, zero(a.half_len))
    s, c = sincos(st.rot)
    rot_mat = SMatrix{2,2}(c, s, -s, c)
    b0 = rot_mat * SVector{2}(-b.half_len, zero(b.half_len)) + st.pos
    b1 = rot_mat * SVector{2}(b.half_len, zero(b.half_len)) + st.pos

    r = b0 - a0
    u = a1 - a0
    v = b1 - b0

    ru = dot(r, u)
    rv = dot(r, v)
    uu = dot(u, u)
    uv = dot(u, v)
    vv = dot(v, v)

    det = uu * vv - uv * uv
    cond = det < eps * uu * vv
    s = IfElse.ifelse(cond, clamp01(ru / uu), clamp01((ru * vv - rv * uv) / det))
    t = IfElse.ifelse(cond, zero(F), clamp01((ru * uv - rv * uu) / det))

    S = clamp01((t*uv + ru)/uu)
    T = clamp01((s*uv - rv)/vv)

    A = a0 + S*u
    B = b0 + T*v

    U = B - A
    dist = sqrt(dot(U, U))
    PA = A + U * a.radius / dist
    PB = B - U * b.radius / dist

    return SVector{5}(
        dist - a.radius - b.radius,
        PA...,
        PB...,
    )
end
