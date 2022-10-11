# Collide

A simple, performant and accurate physics engine implemented in pure Julia using ModelingToolkit.jl and DifferentialEquations.jl. Depends on [PrimitiveCollisions.jl](https://github.com/AayushSabharwal/PrimitiveCollisions.jl) for collisions. To use:

```
] add https://github.com/AayushSabharwal/PrimitiveCollisions.jl
] add https://github.com/AayushSabharwal/Collide.jl
```

A simple world with 2 rectangular entities, and no gravity:

```julia
using Collide, PrimitiveCollisions, StaticArrays

# Create two entities
e = Collide.Entity(
    name = :a,
    shape = PrimitiveCollisions.Rect(1.0, 1.0),
    linear_drag = 0.5
)
e2 = Collide.Entity(
    name = :b,
    shape = PrimitiveCollisions.Rect(1.0, 1.0),
    position = SVector{2}(2.1, -1.2),
    velocity = SVector{2}(-1.0, 0.0),
)

# Add them to the world
world = Collide.World(:w)
push!(world, e)
push!(world, e2)

sim = Collide.Simulation(w, [0., 0.]); # Add this semicolon or face pages of printed content

# Step 10.0 seconds, stopping after exactly 10.0s
Collide.step!(sim, 10.0, true)
```

# Gotchas

- Entities that are intersecting when the collision starts will be stuck inside each other (they "collide from the inside").
- Entities that touch exactly (two circles at a distance equal to the sum of their radii) are also considered to intersect. To simulate things like Newton's cradle, separate the circles by a tiny distance (`eps()`).
- The internals are not yet documented, but not too complicated either.

# Plotting

A simple plotting function, and one to generate a bunch of random bodies inside a box, are available in [this gist](https://gist.github.com/AayushSabharwal/9c76640b960c6726072d42f7ea577e08). Try the following:

```julia
sim, w = manybody(10, 10)
animate(sim, w; tstep = 0.01, frames = 3000)
```
