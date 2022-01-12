# Helper functions for spherical geometry

using Rotations

"""
    c2s(x, y, z) -> [θ, ϕ]
    c2s(xyz::AbstractVector) -> [θ, ϕ]
    c2s(xyz::AbstractMatrix) -> [θs; ϕs]

Convert Cartesian vector `[x, y, z]` (or `xyz`) to spherical unit vector `[θ,
ϕ]`.  For AbstractMatrix `xyz`, each column is converted and the results are
reassmebled into a Matrix using `hcat`.
"""
function c2s(x,y,z)
    θ=atan(y,x)
    ϕ=atan(z, hypot(x,y))
    [θ, ϕ]
end

function c2s(xyz::AbstractVector{<:Real})
    c2s(xyz[1], xyz[2], xyz[3])
end

function c2s(xyz::AbstractMatrix{<:Real})
    @views reduce(hcat, c2s.(xyz[1,:], xyz[2,:], xyz[3,:]))
end

"""
    s2c(θ, ϕ) -> [x, y, z)
    s2c(θϕ::AbstractVector) -> [x, y, z]
    s2c(θϕ::AbstractMaxtrix) -> [xs; ys; zs]

Convert spherical unit vector `[θ, ϕ]` (or `θϕ`) to Cartesian unit vector `[x,
y, z]`.  For AbstractMatrix `θϕ`, each column is converted and the results are
reassmebled into a Matrix using `hcat`.
"""
function s2c(θ,ϕ)
    z, r = sincos(ϕ)
    y, x = sincos(θ) .* r
    [x, y, z]
end

function s2c(θϕ::AbstractVector{<:Real})
    s2c(θϕ[1], θϕ[2])
end

function s2c(θϕ::AbstractMatrix{<:Real})
    @views reduce(hcat, s2c.(θϕ[1,:], θϕ[2,:]))
end

"""
    beamrings(θ, ϕ; nrings=4, dϕ=deg2rad(10/3600)) -> Matrix

Compute the centers of beams that are arranged in concentric rings around the
spherical coordinates `(θ, ϕ)`, both in radians with `θ` being the longitudinal
coordinate and `ϕ` being the latitudinal coordinate.  The number of rings for
which beam centers are returned is given by `nrings`, which defaults to 4.  The
rings are radially separated by `dϕ` radians, which defaults to the equivalent
of 10 arcseconds.  The beams returned do not form a packed hexagon because that
is not possible on a sperical sky, but the number of beams per ring does follow
the hexagonal packing convention of `6n` beams for ring `n` (when `n>0`).  The
so-called boresight beam center at `(θ, ϕ)` is included as well, so 61 beam
centers are returned for `nrings==4`.  More generally,
`nbeams=6*sum(1:nrings)+1` beam centers are returned in a `2xnbeams` Matrix,
where the first row contains the `θ_beam` coordinates and the second row
contains the `ϕ_beam` coordinates.
"""
function beamrings(θ, ϕ; nrings=4, dϕ=deg2rad(10/3600))
    nbeams = 6 * sum(1:nrings) + 1
    θs = Vector{Float64}(undef, nbeams)
    ϕs = Vector{Float64}(undef, nbeams)
    θs[1] = 0.0
    ϕs[1] = π/2
    b=2
    for r in 1:nrings
        nθ = 6r
        dθ = 2π/nθ
        ϕr = π/2 - r * dϕ
        for i in 0:nθ-1
            θs[b] = i*dθ
            ϕs[b] = ϕr
            b += 1
        end
    end
    reduce(hcat, c2s.(Ref(RotZY(θ, π/2-ϕ)) .* s2c.(θs, ϕs)))
end
