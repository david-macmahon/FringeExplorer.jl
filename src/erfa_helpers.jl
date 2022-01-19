# Helper functions to simplify ERFA usage

using ERFA: DAS2R, DJM0, apco13, atciq, atco13, atioq
using RadioInterferometry: dms2rad
using EarthOrientation: getΔUT1, getxp, getyp

const MJD0 = DJM0

"""
    radec2hadec(α, δ, jdutc;
                longitude, latitude, altitude,
                jdutc2=0,
                dut1=getΔUT1(jdutc+jdutc2),
                xp=getxp(jdutc+jdutc2)*DAS2R,
                yp=getxp(jdutc+jdutc2)*DAS2R
               ) -> (hob, δob)

    radec2hadec(αδ::AbstractMatrix, jdutc; kwargs...) -> Matrix[hobs; δobs]

Convert ICRS right ascension `α` and declination `δ` to observed hour angle
`hob` and observed declination `δob` as seen from the location specified by
`longitude`, `latitude`, and `altitude` at the Julian date `jdutc+jdutc2`.
Earth orientation parameters `dut1`, `xp`, and `yp` may be specified if desired,
otherwise values from EarthOrientation for `jdutc+jdutc2` will be used.

If `α` and `δ` are Vectors and `jdutc+jdutc2` is a scalar, the
position-independent astrometric parameters for `jdutc+jdutc2` will be computed
only once and `hob` and `δob` will also be Vectors.  Multiple right ascensions
and declinations given as a 2xN matrix `αδ` consisting of [`α`, `δ`] columns
will also compute the position-independent astrometric parameters only once for
all positions.

The `longitude` and `latitude` parameters should be given in degrees.  They may
be numeric values or Strings in decimal ("dd.ddd") or sexagesimal
("dd:mm:ss.sss") format.

!!! note
    Be sure to use the right units for the various parameters!

| Parameter | Description            | Units   |
|-----------|:-----------------------|:--------|
|         α | Right ascention (ICRS) | radians |
|         δ | Declination (ICRS)     | radians |
|     jdutc | Julian Date (UTC)      | days    |
| longitude | Longitude              | degrees |
|  latitude | Geodetic latitude      | degrees |
|  altitude | Altitude               | meters  |
|    jdutc2 | Julian Date (UTC)      | days    |
|      dut1 | UTC-UT1                | seconds |
|        xp | Polar motion X         | radians |
|        yp | Polar motion Y         | radians |
"""
function radec2hadec(α, δ, jdutc;
                     longitude, latitude, altitude,
                     jdutc2=0,
                     dut1=getΔUT1(jdutc+jdutc2),
                     xp=getxp(jdutc+jdutc2)*DAS2R,
                     yp=getxp(jdutc+jdutc2)*DAS2R,
                     _ignored_kwargs...
                    )
    _aob, _zob, hob, δob, _αob, _eo = atco13(α, δ,
        0, 0, 0, 0,
        jdutc, jdutc2, dut1,
        dms2rad(longitude), dms2rad(latitude), altitude,
        xp, yp,
        0, 0, 0, 0
    )

    (hob, δob)
end

function radec2hadec(α::AbstractVector, δ::AbstractVector, jdutc::Real;
                     longitude, latitude, altitude,
                     jdutc2=0,
                     dut1=getΔUT1(jdutc+jdutc2),
                     xp=getxp(jdutc+jdutc2)*DAS2R,
                     yp=getxp(jdutc+jdutc2)*DAS2R,
                     _ignored_kwargs...
                    )

    # Star-independent astrometry parameters
    astrom, _eo = apco13(jdutc, jdutc2, dut1,
                         dms2rad(longitude), dms2rad(latitude), altitude,
                         xp, yp,
                         0, 0, 0, 0
                        )

    # Transform ICRS to CIRS
    αδis = atciq.(α, δ, 0, 0, 0, 0, Ref(astrom))

    αis = getindex.(αδis, 1)
    δis = getindex.(αδis, 2)

    # Transform CIRS to observed
    azhδαobs = atioq.(αis, δis, Ref(astrom))

    hobs = getindex.(azhδαobs, 3)
    δobs = getindex.(azhδαobs, 4)
    (hobs, δobs)
end

function radec2hadec(αδ::AbstractMatrix, jdutc::Real;
                     longitude, latitude, altitude,
                     jdutc2=0,
                     dut1=getΔUT1(jdutc+jdutc2),
                     xp=getxp(jdutc+jdutc2)*DAS2R,
                     yp=getxp(jdutc+jdutc2)*DAS2R,
                     _ignored_kwargs...
                    )
    @assert size(αδ, 1) == 2 "expected 2xN Matrix"

    hobs, δobs = radec2hadec(αδ[1,:], αδ[2,:], jdutc;
                             longitude, latitude, altitude,
                             jdutc2, dut1, xp, yp
                            )

    [permutedims(hobs); permutedims(δobs)]
end
