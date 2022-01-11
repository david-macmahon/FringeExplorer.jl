module FringeExplorer

using RadioInterferometry
using Symbolics

export ha2t
export t2ha
export td2wdw
export td2wdw!
export hd2wdw
export hd2wdw!

export telinfo
export antpos
include("telinfo.jl")

export s2c
export c2s
export beamrings
include("spherical.jl")

export MJD0
export radec2hadec
include("erfa_helpers.jl")

"""
Rotaional rate of Earth in radians per second
"""
const Ωe = 7.29211538e-5

"""
    ha2t(ha::Real, Ω=Ωe) -> t

Convert hour angle `ha` in radians to time `t` in seconds. `ha == t == 0` at
transit.  `Ω`, the rate of rotation in radians per second, defaults to `Ωe`, the
rotational rate of Earth.
"""
ha2t(ha::Real, Ω=Ωe) = ha / Ω

"""
    t2ha(t::Real, Ω=Ωe) -> ha

Convert time `t` in seconds to hour angle `ha` in radians. `ha == t == 0` at
transit.  `Ω`, the rate of rotation in radians per second, defaults to `Ωe`, the
rotational rate of Earth.
"""
t2ha(t::Real, Ω=Ωe) = t * Ω

# Symbolic variables
const Ω, t, δ = @variables Ω t δ
const Dt = Differential(t)
const h = Ω * t

# Passing zero for longitude requires rotating XYZ position vectors around Z
# such that positive X intersects the local meridian.
const xyzR3uvw = xyz2uvw(h, δ, 0)

# Substitute value to `Ω`
const xyzR3uvw_e = substitute.(xyzR3uvw, [Ω=>Ωe])

# Delay is w component (use 3:3 to keep as a one row Matrix)
const w = xyzR3uvw_e[3:3,:]
# Delay rate is derivative of w component
const dw = expand_derivatives.(Dt.(w))
# Make 2x3 matrix from w and dw
const wdw_m = [w; dw]

# Create expressions for functions that return wdw for given t and δ
const wdw_expr = build_function(wdw_m, t, δ)

"""
    td2wdw(t, δ) -> Matrix
    tδ2wdw(t, δ) -> Matrix

Returns a 2x3 Matrix that can left multiply a 3 element vector (or a 3 row
matrix) of antenna position(s) to get delay `w` in the first element (or row)
and the delay rate `dw` in the second for a source observed at time `t` and
declination `δ`.  The time `t` should be in seconds, with `t<0` before the
source transits, `t==0` when the source transits, and `t>0` after the source
transits.

# Example
```julia
wdw = td2wdw(0, 0) * antpos
```
"""
td2wdw = eval(wdw_expr[1])
const tδ2wdw = td2wdw

"""
    td2wdw!(dst, t, δ)
    tδ2wdw!(dst, t, δ)

Stores a 2x3 Matrix in `dst` (which should already be a 2x3 Matrix) that can
left multiply a 3 element vector (or a 3 row matrix) of antenna position(s) to
get delay `w` in the first element (or row) and the delay rate `dw` in the
second for a source observed at time `t` and declination `δ`.  The time `t`
should be in seconds, with `t<0` before the source transits, `t==0` when the
source transits, and `t>0` after the source transits.
"""
td2wdw! = eval(wdw_expr[2])
const tδ2wdw! = td2wdw!

"""
    hd2wdw(h, δ, Ω=Ωe) -> Matrix
    hδ2wdw(h, δ, Ω=Ωe) -> Matrix

Returns a 2x3 Matrix that can left multiply a 3 element vector (or a 3 row
matrix) of antenna position(s) to get delay `w` in the first element (or row)
and the delay rate `dw` in the second for a source observed at hour angle `h` and
declination `δ`.  The hour angle `h` should be in radians, with `h<0` before the
source transits, `h==0` when the source transits, and `h>0` after the source
transits.  `Ω` is used to convert `h` to time (see [`ha2t`](@ref)).

# Example
```julia
wdw = hd2wdw(ha, dec) * antpos
```
"""
hd2wdw(h, δ, Ω=Ωe) = tdwdw(ha2t(h, Ω), δ)
hδ2wdw = hd2wdw

"""
    hd2wdw!(dst, t, δ, Ω=Ωe)
    hδ2wdw!(dst, t, δ, Ω=Ωe)

Stores a 2x3 Matrix in `dst` (which should already be a 2x3 Matrix) that can
left multiply a 3 element vector (or a 3 row matrix) of antenna position(s) to
get delay `w` in the first element (or row) and the delay rate `dw` in the
second for a source observed at time `t` and declination `δ`.  The time `t`
should be in seconds, with `t<0` before the source transits, `t==0` when the
source transits, and `t>0` after the source transits.  `Ω` is used to convert
`h` to time (see [`ha2t`](@ref)).
"""
hd2wdw!(dst, h, δ, Ω=Ωe) = tdwdw!(dst, ha2t(h, Ω), δ)
hδ2wdw! = hd2wdw!

end # module FringeExplorer
