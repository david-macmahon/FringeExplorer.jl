### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 48c789c8-214e-11eb-0c42-c90206ac9147
# Boilerplate code to pull in required packages in a portable way
begin
	import Pkg
	Pkg.activate(Base.current_project())
	using Plots
	using ERFA
	using LinearAlgebra # for dot product

md"# Phased Array Beam Patterns

This notebook develops a beam pattern analysis for phased arrays of homogenous radio telescopes.
"
end

# ╔═╡ fa66394c-21a1-11eb-1991-fb0fb1181350
md"""## Introduction"""

# ╔═╡ dc0d7392-2175-11eb-0964-5b252639b5fd
md"""For steerable radio telescope arrays, the array's _control system_ directs the antennas' drive systems to orient the antennas' parabolic reflectors such that the center of each of antenna's primary beam is aligned with the location specified by the observer's instructions. This location is referred to as the _pointing center_. Typically, the pointing center is continuously updated to _track_ the apparent location of a celestial source (i.e. object or location) as it moves across the sky.

Most steerable radio telescopes have _altitude-azimuth_ (or _Alt-Az_) drive systems that are steerable in azimuth and elevation (aka altitude).  Some have equatorial (aka polar) drive systems that are steerable in hour angle and declination.  In either case, the phased array beam analysis is similar.

Non-steerable radio telescopes rely on the apparent motion of sources across the sky (e.g. due to the rotation of Earth) to allow observation of multiple pointing centers.  Many of the phased array concepts presented here are relevant for non-steerable radio telescopes as well, but non-steerable radio telescopes often have non-circular primary beam patterns which can complicate such analysis.

The analysis presented here is for a radio telescope array whose antennas have a circular primary beam, though the concepts should apply to non-circular primary beams as well.  This analysis also assumes that the radio telescopes are steerable in azimuth and elevation, but the end result is basically the same regardless of the type of drive system.
"""

# ╔═╡ b2c21684-21a5-11eb-0d02-853dfc69164c
md"## Phased Array Theory of Operation"

# ╔═╡ ceed38a2-231f-11eb-1f9a-cf730d080095
md"""The voltage produced by each radio telescope is proportional to the sum of the radio frequency energies emitted by the all the sources within the primary beam (and sidelobes) plus the uncorrelated system noise of the antenna and signal path."""

# ╔═╡ 7f734cbc-219b-11eb-1471-594ab489e084
md"When operating in _phased array_ mode, the voltages from the antennas are aligned in both time and phase relative to a desired sky position and the array's reference position. This desired sky position is known as the _phase center_. Although not strictly required, the phase center is almost always located within the primary beam. When the phase center happens to be the same as the pointing center (i.e. the center of the primary beam), the array is said to be _phased to boresight_.
"

# ╔═╡ fd656238-2176-11eb-3a4e-7dbd106b836b
md"""The time and phase aligned signals from all the antennas are then added together to produce a _phased array sum_ (aka _coherent beam_).  This results in _constructive interference_ at the phase center where the correlated energy (i.e. the energy from celestial sources, not the energy from system noise) from the antennas add together in-phase for maximum gain while the uncorrelated energy (aka system noise) adds together randomly and the correlated energy from elsewhere in the primary beam adds together with sub-optimal phasing that varies across the primary beam. Determining how much gain each primary beam location experiences due to summing is the main objective of the phased array beam pattern analysis.

The time alignment depends on the _geometric delays_ from the antennas to the array reference location as the RF wavefront propagates across the array. These geometric delays differ depending on phase center location and array geometry. The signal from each antenna also experiences a finite _fixed delay_ that generally differs between antennas (e.g. due to cable length variations). The geometric delays vary as the phase center moves relative to the array whereas the fixed delays, as their name implies, do not generally vary.  Geometric delays can be calculated when the phase center and array geometry are known, but fixed delays must be measured in situ or given a priori.  Other delays, such as those induced by atmospheric variations over the array, can further complicate things but this analysis will ignore these other sources of delay and delay variation.

Once the geometric and fixed delays are known and compensated for, the signals must still be phase aligned before they will add constructively.  The signal from an antenna can be delay corrected for a single delay value only.  This delay value is the geometric delay for the _delay center_ (generally the same location as the phase center) plus the antenna's fixed delay. Delay corrections are generally performed in discrete steps, which means that there is usually a small residual delay that remains uncorrected.  Phase is the product of frequency and time, ``\phi = \omega \tau``, so the residual delay (a time offset) results in a frequency dependent phase offset.  This phase offset can be removed on a per channel basis after the signal has been channelized into multiple frequency components (e.g. by an FFT or PFB). Locations other than the delay center also have their own residual delays that vary depending on their relationship to the delay center.  Additionally, each antenna has a unique instrumental phase that results from many different factors and must be measured in situ.  Determination of the instrumental phases is beyond the scope of this analysis, so we assume here that the instrumental phases are zero or have already been corrected for.
"""

# ╔═╡ e7853568-21a5-11eb-13de-f17b68612831
md"""Phasing the antennas to a specific phase center requires determination of the geometric delays.  The concept for this is relatively straightforward.  The geometric delay is proportional to the geometric distance that the planar wavefront travels from antenna to array reference location. This distance is the projection of the baseline vector ``\vec{b}``, defined to be the vector from the array reference position to the antenna position, onto the unit vector ``\hat{s}``, which points in the direction of the delay center from the array reference position.  This can be computed by taking the dot product of ``\vec{b}`` and ``\hat{s}``, which is written as ``\vec{b}\cdot\hat{s}``.  The array reference position is usually treated as the origin of the array reference frame so the baseline vector ``\vec{b}`` is the same as the antenna position in the array reference frame.  To take the dot product, both ``\vec{b}`` and ``\hat{s}`` must be represented in the same frame."""

# ╔═╡ 190282de-21a2-11eb-176e-178fd029cb6d
md"""## Primary Beam Coordinates"""

# ╔═╡ c2295c46-219f-11eb-2d92-cfdeba1d942d
md"""The analysis starts by exploring the coordinates of various points within the primary beam of the telescope. The primary beam is modeled as a circle with an angular diameter equal to the full width half maximum (FWHM) of a single radio telescope. For convenience, this circle is inscribed in a square region of sky which is represented by a two dimensional grid of points with the central point corresponding to the center of the primary beam, i.e. the pointing center. To keep things simple (or at least less complicated), this analysis uses a _flat sky_ approximation.
"""

# ╔═╡ 54a02260-219d-11eb-17b6-17b0e0bc9d8d
md"The axes of this grid correspond to cross-elevation (horizontal, aka `xel`) and elevation (vertical, aka `el`) as seen from the telescope's perspective. The origin of the grid is the center of the primary beam. The unit steps in the cross-elevation and elevations directions are `dxel` and `del`, respectively, though they are sometimes shortened here to `dx` and `de`. The coordinate frame of this grid is referred to here as the _primary beam frame_. A primary beam with FWHM of 2 (arbitrary angular units) and `dxel`, `del` values of 0.25 (arbitrary angular units) can depicted as: 
"

# ╔═╡ bc7cb314-2150-11eb-2900-9b1c7ab0a8d9
begin
	plot(cos, sin, 0:0.01:2π, aspect_ratio=:equal,
		legend=false, xflip=true,
		xlims=(-1.1,1.1), ylims=(-1.1,1.1)
	)
	scatter!(vec([(x,y) for y=-1:0.25:1, x=-1:0.25:1]), ms=2)
	title!("Primary Beam")
	xlabel!("Cross Elevation (dxel=0.25)")
	ylabel!("Relative Elevation (del=0.25)")
end

# ╔═╡ 2f22727e-22ef-11eb-379f-fd85d1a81d55
md"""For reasons that will be clear later, we have opted to orient the cross elevation (horizontal) axis with positive values to the left of the origin and negative values to the right."""

# ╔═╡ a018417a-2227-11eb-2f6d-b7460d2e140d
md"""After compensating for geometric delay at the delay center, the residual delay at the delay center is zero.  The rest of the antenna's primary beam will have a delay gradient that depends on the geometry of the delay center relative to the antenna location and array reference location.  To compute the delay for a given point, we must create a unit vector, ``\hat{s}``, in the direction of the point and take the dot product of that and the baseline vector ``\vec{b}``.  We are free to choose the coordinate frame in which ``\hat{s}`` is represented.  Usually that frame will be different from the frame in which the antenna positions are expressed, so before taking the dot product we need to ensure we convert the antenna positions into the frame of ``\hat{s}`` or vice versa.  For beam pattern analysis we will typically have more ``\hat{s}`` vectors than antennas, so it will be more efficient to transform the antenna positions into the frame of the ``\hat{s}`` vectors."""

# ╔═╡ 9eb14aca-2323-11eb-3543-33f08ff00e2c
md"""As described above, the coordinate frame for ``\hat{s}`` is most naturally expressed in terms of spherical coordinates of ``(\theta, \phi, 1)``, where ``\theta`` is ``\sphericalangle xel``, ``\phi`` is ``\sphericalangle el``, and ``1`` is the unit length of ``\hat{s}``. The ``\hat{s}`` vector needs to be expressed in a rectilinear form so that we can take the dot product with rectilinear baseline vector ``\vec{b}``. The rectalinear form of ``\hat{s}`` is also known as _direction cosines_. The _Essential Routines for Fundamental Astronomy_ (ERFA) package provides a spherical to Cartesian function, ``\mathtt{s2c}(\theta,\phi)``, that we can use to calculate ``\hat{s}`` for any values of `xel` and `el`.  Here is how to calculate ``\hat{s}`` for the pointing center and two other points (one on each axis):"""

# ╔═╡ 5b49716c-2327-11eb-2df4-37ceb3707973
begin
	dxel = 0.25 * ERFA.DD2R
	del  = 0.25 * ERFA.DD2R
	[ERFA.s2c(θ, ϕ) for (θ, ϕ) in [(0,0), (dxel,0), (0,del)]]
end

# ╔═╡ a85178f6-2327-11eb-2535-c1f23b6a5484
md"""As you can see, `s2c()` takes two angles (NB: in radians) and returns a rectilinear vector that, in our case, corresponds to `(s, xel, el)`, where `s` is the direction of the source (away from the viewer), `xel` is the cross elevation direction, and `el` is the elevation direction. We choose `xel` to be positive to the left of the origin and `el` positive above the origin to make this a right-handed system because the antenna positions are also specified in a right handed system."""

# ╔═╡ a0fd25d6-2335-11eb-223d-fb1a8be4a076
md"""## Antenna Position Coordinates"""

# ╔═╡ ac2e3c86-222b-11eb-0297-f1aa2b01b60d
md"""Antenna positions are often represented in a topocentric (East, North, Up) or ENU frame, which is a right-handed system.  Antenna positions can also be represented in an _Earth Centered Earth Fixed_ (ECEF) frame known as the _International Terrestrial Reference Frame_ (ITRF).  The ECEF frame can be transformed to the ENU frame used in the analysis here, so we will not treat it separately.  Antenna positions are usually given in meters, but other units are sometimes used as well (e.g. nanoseconds or wavelengths at a given frequency).  For our purposes, we will ultimately need to convert the units to time and/or a frequency dependent phase."""

# ╔═╡ 5493615c-2336-11eb-2c36-0772eadafb39
md"""## Coordinate Frame Transformations"""

# ╔═╡ 396e99c8-2336-11eb-24a9-fddf80b160b5
md"""The (s, xel, el) primary beam frame described above is beam-centric in that it does not depend on the location of the pointing center on the sky.  On the other hand, the transformation from the topocentric ENU frame to the beam-centric (s, xel, el) frame does depend on the location of the pointing center."""

# ╔═╡ ff98a0d8-2252-11eb-1056-59c804698c62
md"""The transformation from the ENU frame to the (s, xel, el) frame can be represented as a series of two rotations around various axes of the frame itself. Transforming from the ENU frame to the (s, xel, el) frame can be performed by:

  1. Rotating the (E,N,U) frame around the U (Up) axis in the clockwise direction (as viewed looking down to the origin) by the topocentric azimuth angle of the pointing center minus 90 degrees.  This results in an intermediate (E',N',U) frame.
  2. Rotating the intermediate (E',N',U) frame around the N' axis in the clockwise direction (as viewed looking at the origin from positive N') by the elevation angle.  This results in a (E", N', U') frame that is equivalent to the (s, xel, el) frame in which our ``\hat{s}`` vectors are expressed.
"""

# ╔═╡ f4f6ab3a-222c-11eb-0564-17ceb3dde858
md"""Each of these rotations can be represented by a 3x3 rotation matrix. These matrices can be multiplied together, in the proper order, to produce a single composite 3x3 rotation matrix that transforms the ENU frame into the beam-centric (s, xel, el) frame for the azimuth and elevation of a given the pointing center.  The ENU antenna position(s) represented as a 3 element column vector (or 3xM matrix of M antenna positions) can be left multiplied by this composite rotation matrix to produce the antenna position(s) in the (s, xel, el) frame.  Here we define an `enu2beam(az, el)` function that returns this composite rotation matrix for a pointing center at azimuth `az` and elevation `el` (both given in degrees):"""

# ╔═╡ 57172622-2325-11eb-21f5-87aba95fac55
begin
	"""
	Return a 3x3 rotation matrix for ENU to beam transformation for the
	pointing center at azimuth `az` and elevation `el` (both given in degrees).
	"""
	function enu2beam(az, el)
		ERFA.ry(-el*ERFA.DD2R,
			ERFA.rz(π/2 - az*ERFA.DD2R,
				[1.0 0.0 0.0
				 0.0 1.0 0.0
				 0.0 0.0 1.0]
			)
		)
	end

	# Version that takes a Tuple and auto-splats it
	enu2beam(azel::Tuple{Real, Real}) = enu2beam(azel...)
end

# ╔═╡ 610ca96c-2330-11eb-1869-593e53e74f42
md"""Let's try it out on some basic test cases!"""

# ╔═╡ 6096037a-2330-11eb-1e75-c360e0e0c494
enu2beam(0, 0) # north horizon

# ╔═╡ 76617054-2332-11eb-15f1-91a9f5207081
enu2beam(90, 0) # east horizon

# ╔═╡ 440f67c2-2333-11eb-0975-f362da0edde5
enu2beam(0, 90) # zenith from north

# ╔═╡ b26834ae-2333-11eb-13da-7d924b886cec
enu2beam(180, 90) # zenith from south

# ╔═╡ e902f9fe-2333-11eb-0bb5-9949d2f05225
md"""Notice that in the last example (zenith from south), the sign of `xel` and `el` will be opposite from the example that precedes it (zenith from north) even though both have the local zenith as the pointing center."""

# ╔═╡ e8db8432-2333-11eb-3195-877803e736ba
md"""## Putting It All Together"""

# ╔═╡ e8b55f50-2333-11eb-0f6c-abc07cdfd3cc
md"""We now have the tools to start building and putting together all the various parts of the analysis:
  - Grid of points spanning the angular region of interest
  - Unit ``\hat{s}`` vectors for each of those points
  - Antenna positions in ENU form
  - Conversion from ENU frame to beam frame for arbitrary pointing center
  - Dot product to compute distance/delay from antenna to array reference position
"""

# ╔═╡ d620f9aa-23ba-11eb-11a2-75d99569a381
md"""Here we create a grid of points spanning the angular region of interest and unit ``\hat{s}`` vectors corresponding to those points:"""

# ╔═╡ 165f5502-23b6-11eb-3824-dde609547349
npoints = 65

# ╔═╡ 1f455e96-2544-11eb-3f08-4b013d6b71f1
span_degrees = 4.0

# ╔═╡ 11735a72-23b7-11eb-0c75-39f94e10ba04
xels = range(-span_degrees/2, +span_degrees/2, length=npoints)

# ╔═╡ 4697cad0-23b7-11eb-17fd-d73044fab8bb
els = range(-span_degrees/2, +span_degrees/2, length=npoints)

# ╔═╡ 7ba15bf6-23b7-11eb-3355-9932eae00a13
grid_coords = collect(Iterators.product(xels, els))

# ╔═╡ a121be52-23b7-11eb-3fa0-f7b7cef14f66
# We could also define an `ERFA.s2c(t::Tuple)` method and
# use broadcasting, but using `map()` is good enough for now.
grid_svecs = map(p->ERFA.s2c((p.*ERFA.DD2R)...), grid_coords)

# ╔═╡ babfe764-23ba-11eb-02a2-714e5c641733
md"""We start with one antenna which you are free to move around (up to 1 km in each direction):

$(@bind user_east html"<input type=range min=-100000 max=100000 value=3000 style='width: 75%'> East</br>")
$(@bind user_north html"<input type=range min=-100000 max=100000 value=0 style='width: 75%'> North</br>")
$(@bind user_up html"<input type=range min=-100000 max=100000 value=0 style='width: 75%'> Up")
"""

# ╔═╡ bf0e1a08-23e1-11eb-22e3-2f9770941b9f
user_antpos_enu = [user_east, user_north, user_up] ./ 100

# ╔═╡ 34d02b10-23e1-11eb-1ef6-3540d2475004
md"""And of course you get to point the antenna too:

$(@bind user_az html"<input type=range min=0 max=359 value=90 style='width: 75%'> Azimuth</br>")
$(@bind user_el html"<input type=range min=0 max=90 value=0 style='width: 75%'> Elevation")
"""

# ╔═╡ ba454464-23ba-11eb-08f7-bb9bd85194e8
user_pointing_azel = (user_az, user_el)

# ╔═╡ 62a0f894-23d7-11eb-347f-c3c1de1f41a8
md"""Now we can transform your specified (E, N, U) antenna position from the ENU frame to the (s, xel, el) frame of your specified (az, el) pointing center."""

# ╔═╡ 18d78b34-23d6-11eb-33f1-8df5a4066d65
user_enu2beam = enu2beam(user_pointing_azel)

# ╔═╡ e899c044-23d6-11eb-2780-a1fe38bec119
user_antpos_pointing = user_enu2beam * user_antpos_enu

# ╔═╡ ece4a6b8-23c3-11eb-38e7-bd248d0f0bc5
md"""With the antenna position in the same frame as the unit ``\hat{s}`` vectors, we can take the dot products to calculate the distances from the antenna to the array reference position (i.e. the origin of the ENU frame) for each point in the primary beam.  We are more interested in the relative distances across the beam rather than the absolute distances.  To calculate the relative distances we subtract each unit vector ``\hat{s}`` from unit vector ``\hat{d}`` that points in the direction of the delay center (assumed here to be the same as the pointing center).  These relative distances can be plotted as a contour plot.
"""

# ╔═╡ a2bf4b28-23d5-11eb-1ef1-1170b39fcb9e
reldists = map(s->dot(ERFA.s2c(0,0)-s, user_antpos_pointing), grid_svecs);

# ╔═╡ 601ffcd6-23c4-11eb-2012-01bb72855b38
begin
heatmap(xels, els, reldists', xflip=true,
	aspect_ratio=:equal, fillcolor=:lightrainbow,
	xlims=extrema(xels), ylims=extrema(els),
	title="Antenna to Origin: Distance\nrelative to (az=$(user_az), el=$(user_el))",
	xlabel="Cross Elevation (degrees)",
	ylabel="Elevation (degrees)"
)
contour!(xels, els, reldists', xflip=true,
	aspect_ratio=:equal, fill=false, clabels=true,
	xlims=extrema(xels), ylims=extrema(els),
	fillcolor=:lightrainbow, linecolor=:black,
	levels=0.01:0.01:1,
	colorbar_entry=false)
end

# ╔═╡ a37ac5d2-2483-11eb-3761-c18b249d25c1
md"""For a given frequency, these relative distances can be converted to relative phases which can also be plotted as a heatmap plot:

$(@bind user_mhz html"<input type=range min=500 max=2000 value=1420 style='width: 75%'> Frequency (MHz)")
"""

# ╔═╡ af50e5c2-2487-11eb-1da9-add2857380a8
md"""Frequency: $(user_mhz) MHz"""

# ╔═╡ 13c703ec-2488-11eb-06e2-a93d6fc4ada1
begin
	# ERFA.CMPS is speed of light in meters per second
	λ = ERFA.CMPS / user_mhz / 1e6
	relphases = rem2pi.(2π .* reldists ./ λ, RoundNearest) .* ERFA.DR2D
end;

# ╔═╡ de5be0d6-2489-11eb-21b1-e3218333e30d
begin
	plt=heatmap(xels, els, relphases', xflip=true,
		aspect_ratio=:equal, color=:hsv,
		xlims=extrema(xels), ylims=extrema(els), clims=(-180,180), levels=71,
		title="Antenna to Origin: Degrees of Phase\nrelative to az=$(user_az), el=$(user_el) at $(user_mhz) MHz (λ=$(round(1e5λ)/1e3) cm)\n",
		xlabel="Cross Elevation (degrees)",
		ylabel="Elevation (degrees)"
	)
	# Plot 5 degree contours if phases are within +/- 90 degrees
	# (contouring around phase wraps gets ugly)
	if maximum(abs.(relphases)) <= 90
		contour!(xels, els, relphases', xflip=true,
			aspect_ratio=:equal, fill=false, clabels=true,
			levels=-90:5:90, linecolor=:black,
			colorbar_entry=false
		)
	end
	plt
end

# ╔═╡ 7324d316-24ae-11eb-2112-a51606afab74
if maximum(abs.(relphases)) <= 90
	md"""(Phases within ±90 degrees: Showing 5 degree contour lines.)"""
else
	md"""(Phases greater than ±90 degrees: Contour lines are not being shown.)"""
end

# ╔═╡ 22338236-24b3-11eb-02aa-0bc76375dfd3
md"""## Adding Multiple Antennas"""

# ╔═╡ 4196dd3a-24b3-11eb-2e0d-650a39dc4bb6
md"""Now that we can compute phases across the primary beam of an antenna relative to the array center for a given pointing/delay center at a given frequency, the next thing to do is add more antennas to the array and then form a coherent beam by summing them.  With multiple antennas, we can delay correct them for a common delay/pointing center, calculate their phases across the primary beam, convert these phases to unit vectors in the complex plane (i.e. complex numbers), and add corresponding complex numbers together to create the complex beam pattern of the coherent beam created by adding all the delay corrected signals together."""

# ╔═╡ aef63f92-24b3-11eb-1a6c-79b47130a4a9
md"""We will use four new antennas, a new pointing center, and new frequency.  To keep the UI simple, the new antennas will have an `Up` coordinate of zero, but you can still move them East/West and North/South."""

# ╔═╡ f71a7bb4-26f5-11eb-2c17-cb548d2879c5
@bind ants_reset html"<input type=button value='Reset Antennas'>"

# ╔═╡ 129dca5e-2588-11eb-2e0d-5941534cfde0
begin
	ants_reset
	md"""
	$(@bind ant1_east html"<input type=range min=-100000 max=100000 value=3000 style='width: 75%'> Ant 1 East</br>")
	$(@bind ant1_north html"<input type=range min=-100000 max=100000 value=0 style='width: 75%'> Ant 1 North</br>")
	$(@bind ant2_east html"<input type=range min=-100000 max=100000 value=0 style='width: 75%'> Ant 2 East</br>")
	$(@bind ant2_north html"<input type=range min=-100000 max=100000 value=3000 style='width: 75%'> Ant 2 North</br>")
	$(@bind ant3_east html"<input type=range min=-100000 max=100000 value=-3000 style='width: 75%'> Ant 3 East</br>")
	$(@bind ant3_north html"<input type=range min=-100000 max=100000 value=0 style='width: 75%'> Ant 3 North</br>")
	$(@bind ant4_east html"<input type=range min=-100000 max=100000 value=-0 style='width: 75%'> Ant 4 East</br>")
	$(@bind ant4_north html"<input type=range min=-100000 max=100000 value=-3000 style='width: 75%'> Ant 4 North</br>")
	$(@bind ants_az html"<input type=range min=0 max=359 value=90 style='width: 75%'> Azimuth</br>")
	$(@bind ants_el html"<input type=range min=0 max=90 value=0 style='width: 75%'> Elevation")
	$(@bind ants_mhz html"<input type=range min=500 max=2000 value=1420 style='width: 75%'> Frequency")
	"""
end

# ╔═╡ 15fea64e-24ba-11eb-09cf-4d47980c7fb5
begin
	ants_enu2beam = enu2beam(ants_az, ants_el)

	ants_antpos = hcat(
		[ant1_east, ant1_north, 0],
	    [ant2_east, ant2_north, 0],
	    [ant3_east, ant3_north, 0],
	    [ant4_east, ant4_north, 0]) / 100

	ants_antpos_pointing = ants_enu2beam * ants_antpos

	ants_reldists = map(eachcol(ants_antpos_pointing)) do antpos_pointing
		map(s->dot(ERFA.s2c(0,0)-s, antpos_pointing), grid_svecs)
	end


	ants_λ = ERFA.CMPS / ants_mhz / 1e6
	ants_relphases = map(ants_reldists) do reldists
		rem2pi.(2π .* reldists ./ ants_λ, RoundNearest)
	end

	ant_plots = map(enumerate(ants_relphases)) do (a, relphases)
		heatmap(xels, els, relphases' .* ERFA.DR2D, xflip=true,
			aspect_ratio=:equal, color=:hsv,
			xlims=extrema(xels), ylims=extrema(els), clims=(-180,180), levels=71,
			title="Ant $(a) Phases",
			#xlabel="Cross Elevation (degrees)",
			#ylabel="Elevation (degrees)"
		)
	end

	plot(ant_plots...)
end

# ╔═╡ 4d3b4aa8-251d-11eb-0937-05fa55ee8cbf
begin
	antsz = map(ants_relphases) do relphases
		exp.(1im*relphases)
	end

	# Calculate primary beam response
	local ν = ants_mhz/1e3 # GHz
	local θb = 87.85 / 60 / ν * ERFA.DD2R
	local ρ = map(t->ERFA.seps(0,0,(t.*ERFA.DD2R)...), grid_coords)
	local pb_amp = cos.(1.189π*ρ/θb)./(1 .- 4(1.189ρ/θb).^2)

	beamz = reduce(+, antsz) .* pb_amp
	mag_beamz = abs.(beamz)

	amp_plot = heatmap(xels, els, 10*log10.(mag_beamz'), xflip=true,
		aspect_ratio=:equal, fillcolor=:thermal,
		xlims=extrema(xels), ylims=extrema(els), clims=(-50,12),
		title="Coherent Beam\nAmplitude (dB)\naz=$(ants_az)°, el=$(ants_el)°",
		titlefontsize=12, 
	)
	contour!(xels, els, 10*log10.(mag_beamz'), xflip=true,
		clims=(-50,12), levels=-40:10:40,
		aspect_ratio=:equal, fill=false, clabels=true,
		fillcolor=:thermal, linecolor=:black,
		colorbar_entry=true
	)

	pha_plot = heatmap(xels, els, angle.(beamz)' .* ERFA.DR2D, xflip=true,
		aspect_ratio=:equal, color=:hsv,
		xlims=extrema(xels), ylims=extrema(els), clims=(-180,180), levels=71,
		title="Coherent Beam\nPhase (degrees)\n freq $(ants_mhz) MHz",
		titlefontsize=12
	)

	plot(amp_plot, pha_plot, size=(680, 350))
end

# ╔═╡ 31420f58-2545-11eb-13cc-3365a66ed1f2
begin
	local lim = 1.1 * maximum(abs.(ants_antpos))

	antpos_plt = scatter(map(v->(v[1],v[2]), eachcol(ants_antpos)),
		lims=(-lim,+lim), border=:frame, aspect_ratio=:equal, legend=false,
		color=palette(:default,4)[1:4], title="Antenna Positions",
		series_annotation=map(s->text(s,8), ["1","2","3","4"]),
		ms=10
	)

	local azel_r = cosd(ants_el)
	local azel_x = azel_r * cosd(90 - ants_az)
	local azel_y = azel_r * sind(90 - ants_az)

	azel_plt = plot(cos, sin, 0:0.02:2π, aspect_ratio=:equal,
		ticks=false, border=:none, color=:black,
		grid=false, legend=false, lims=(-1.1,+1.1),
		annotations=[(0,1.1,"N"), (1.1,0,"E"), (0,-1.1,"S"), (-1.1,0,"W")],
		title="Pointing/Delay Center"
	)
	scatter!((azel_x, azel_y), ms=6, color=:black)

	plot(antpos_plt, azel_plt, size=(600,350))
end

# ╔═╡ b9e16b1a-23ba-11eb-08ad-8fc21f98dca9
md"""## Primary Beam

The primary beam pattern incorporated in the coherent beam pattern above is based on the primary beam model of the 13.5 meter MeerKAT antennas:
"""

# ╔═╡ 3e4354d4-2151-11eb-384f-7706d85ce220
begin
	local ν = ants_mhz/1e3 # GHz
	local θb = 87.85 / 60 / ν * ERFA.DD2R
	local ρ = map(t->ERFA.seps(0,0,(t.*ERFA.DD2R)...), grid_coords)
	local pb_amp = (cos.(1.189π*ρ/θb)./(1 .- 4(1.189ρ/θb).^2)).^2
	heatmap(xels, els, 10*log10.(pb_amp),
		xlims=extrema(xels), ylims=extrema(els), clims=(-50,0), levels=4,
		aspect_ratio=:equal, fill=true, clabels=true,
		fillcolor=:thermal, linecolor=:black,
		title="Primary Beam Power Pattern (dBc)\n$(ants_mhz) MHz"
	)
	contour!(xels, els, 10*log10.(pb_amp),
		clims=(-50,0), levels=-40:10:40,
		aspect_ratio=:equal, fill=false, clabels=true,
		fillcolor=:thermal, linecolor=:black,
		colorbar_entry=false
	)
end	

# ╔═╡ 25a3e7f8-26ff-11eb-292c-2bd4bb7858fd
md"""## What's Next?

The coherent phased array beam pattern analysis techniques developed here lay the foundation for real world analyses, but some details have been glossed over or ignored in the interest of simplicity.  The next steps would be to delve into these details.  This final section describes some of these details.

The techniques presented here have assumed that the pointing center, delay center, and phase center are all the same point.  These techniques need to be enhanced to allow these three points to be independently specified.

Another related task would be to derive the antenna-based phase offsets for a phase center at a given celestial position within the primary beam for a specific radio telescope array.  That analysis could be extended to include a determination of the rate of change in the antenna-based phases as the celestial position moves across the sky.

For frequency domain beam forming, a single phase offset is applied to each frequency channel.  Usually the phase offset will be calculated for the center frequency of each channel.  The frequencies at the edge of the channel will have a residual phase offset since the phase offset calculated for the center frequency of the channel will not be the same as the actual phase offset at the edge frequencies.  This residual phase offset (actually a phase slope across the frequency channel) affects phasing efficiency.

When adding the antennas together earlier, we gave each antenna equal significance or weight.  This assumes that the antenna have equal sensitivity and equal system temperatures.  In a real world application that is not the case and the antennas will need to be weighted differently.  There may also be a benefit to weighting antennas differently depending on where in the array they are.
"""

# ╔═╡ Cell order:
# ╟─48c789c8-214e-11eb-0c42-c90206ac9147
# ╟─fa66394c-21a1-11eb-1991-fb0fb1181350
# ╟─dc0d7392-2175-11eb-0964-5b252639b5fd
# ╟─b2c21684-21a5-11eb-0d02-853dfc69164c
# ╟─ceed38a2-231f-11eb-1f9a-cf730d080095
# ╟─7f734cbc-219b-11eb-1471-594ab489e084
# ╟─fd656238-2176-11eb-3a4e-7dbd106b836b
# ╟─e7853568-21a5-11eb-13de-f17b68612831
# ╟─190282de-21a2-11eb-176e-178fd029cb6d
# ╟─c2295c46-219f-11eb-2d92-cfdeba1d942d
# ╟─54a02260-219d-11eb-17b6-17b0e0bc9d8d
# ╟─bc7cb314-2150-11eb-2900-9b1c7ab0a8d9
# ╟─2f22727e-22ef-11eb-379f-fd85d1a81d55
# ╟─a018417a-2227-11eb-2f6d-b7460d2e140d
# ╟─9eb14aca-2323-11eb-3543-33f08ff00e2c
# ╠═5b49716c-2327-11eb-2df4-37ceb3707973
# ╟─a85178f6-2327-11eb-2535-c1f23b6a5484
# ╟─a0fd25d6-2335-11eb-223d-fb1a8be4a076
# ╟─ac2e3c86-222b-11eb-0297-f1aa2b01b60d
# ╟─5493615c-2336-11eb-2c36-0772eadafb39
# ╟─396e99c8-2336-11eb-24a9-fddf80b160b5
# ╟─ff98a0d8-2252-11eb-1056-59c804698c62
# ╟─f4f6ab3a-222c-11eb-0564-17ceb3dde858
# ╠═57172622-2325-11eb-21f5-87aba95fac55
# ╟─610ca96c-2330-11eb-1869-593e53e74f42
# ╠═6096037a-2330-11eb-1e75-c360e0e0c494
# ╠═76617054-2332-11eb-15f1-91a9f5207081
# ╠═440f67c2-2333-11eb-0975-f362da0edde5
# ╠═b26834ae-2333-11eb-13da-7d924b886cec
# ╟─e902f9fe-2333-11eb-0bb5-9949d2f05225
# ╟─e8db8432-2333-11eb-3195-877803e736ba
# ╟─e8b55f50-2333-11eb-0f6c-abc07cdfd3cc
# ╟─d620f9aa-23ba-11eb-11a2-75d99569a381
# ╟─165f5502-23b6-11eb-3824-dde609547349
# ╟─1f455e96-2544-11eb-3f08-4b013d6b71f1
# ╠═11735a72-23b7-11eb-0c75-39f94e10ba04
# ╠═4697cad0-23b7-11eb-17fd-d73044fab8bb
# ╠═7ba15bf6-23b7-11eb-3355-9932eae00a13
# ╠═a121be52-23b7-11eb-3fa0-f7b7cef14f66
# ╟─babfe764-23ba-11eb-02a2-714e5c641733
# ╟─bf0e1a08-23e1-11eb-22e3-2f9770941b9f
# ╟─34d02b10-23e1-11eb-1ef6-3540d2475004
# ╟─ba454464-23ba-11eb-08f7-bb9bd85194e8
# ╟─62a0f894-23d7-11eb-347f-c3c1de1f41a8
# ╠═18d78b34-23d6-11eb-33f1-8df5a4066d65
# ╠═e899c044-23d6-11eb-2780-a1fe38bec119
# ╟─ece4a6b8-23c3-11eb-38e7-bd248d0f0bc5
# ╠═a2bf4b28-23d5-11eb-1ef1-1170b39fcb9e
# ╟─601ffcd6-23c4-11eb-2012-01bb72855b38
# ╟─a37ac5d2-2483-11eb-3761-c18b249d25c1
# ╟─af50e5c2-2487-11eb-1da9-add2857380a8
# ╠═13c703ec-2488-11eb-06e2-a93d6fc4ada1
# ╟─de5be0d6-2489-11eb-21b1-e3218333e30d
# ╟─7324d316-24ae-11eb-2112-a51606afab74
# ╟─22338236-24b3-11eb-02aa-0bc76375dfd3
# ╟─4196dd3a-24b3-11eb-2e0d-650a39dc4bb6
# ╟─aef63f92-24b3-11eb-1a6c-79b47130a4a9
# ╟─15fea64e-24ba-11eb-09cf-4d47980c7fb5
# ╟─129dca5e-2588-11eb-2e0d-5941534cfde0
# ╟─4d3b4aa8-251d-11eb-0937-05fa55ee8cbf
# ╟─31420f58-2545-11eb-13cc-3365a66ed1f2
# ╟─f71a7bb4-26f5-11eb-2c17-cb548d2879c5
# ╟─b9e16b1a-23ba-11eb-08ad-8fc21f98dca9
# ╠═3e4354d4-2151-11eb-384f-7706d85ce220
# ╟─25a3e7f8-26ff-11eb-292c-2bd4bb7858fd
