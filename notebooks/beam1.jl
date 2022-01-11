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
md"""For steerable radio telecscope arrays, the array's _control system_ directs the antennas' drive systems to orient the antennas' parabolic reflectors such that the center of each of antenna's primary beam is aligned with the location specified by the observer's instructions. This location is referred to as the _pointing center_. Typically, the pointing center is continuously updated to _track_ the apparent location of a celestial source (i.e. object or location) as it moves across the sky.

Most steerable radio telescopes have _altitude-azimuth_ (or _Alt-Az_) drive systems that are steerable in azimuth and elevation (aka altitude).  Some have equitorial (aka polar) drive systems that are steerable in hour angle and declination.  In either case, the phased array beam analysis is similar.

Non-steerable radio telescopes rely on the apparent motion of sources across the sky (e.g. due to the rotation of Earth) to allow observation of multiple pointing centers.  Many of the phased array concepts presented here are relevant for non-steerable readio telescopes as well, but non-steerable radio telescopes often have non-circular primary beam patterns which can complicate such analysis.

The anaysis presented here is for a radio telescope array whose antennas have a circular primary beam, though the concepts should apply to non-circular primary beams as well.  This analysis also assumes that the radio telescopes are steerable in azimuth and elevation, but the end result is basically the same regardless of the type of drive system.
"""

# ╔═╡ b2c21684-21a5-11eb-0d02-853dfc69164c
md"## Phased Array Theory of Operation"

# ╔═╡ ceed38a2-231f-11eb-1f9a-cf730d080095
md"""The voltage produced by each radio telescope is proportional to the sum of the radio frequency energies emitted by the all the sources within the primary beam (and sidelobes) plus the uncorrelated system noise of the antenna and signal path."""

# ╔═╡ 7f734cbc-219b-11eb-1471-594ab489e084
md"When operating in _phased array_ mode, the voltages from the antennas are aligned in both time and phase relative to a desired sky position and the array's reference position. This desired sky position is known as the _phase center_. Although not strictly required, the phase center is almost always located within the primary beam. When the phase center happens to be the same as the pointing center (i.e. the center of the primary beam), the array is said to be _phased to boresight_.
"

# ╔═╡ fd656238-2176-11eb-3a4e-7dbd106b836b
md"""The time and phase aligned signals from all the antennas are then added together to produce a _phased array sum_ (aka _coherent beam_).  This results in _constructive interference_ at the phase center where the correlated energy (i.e. the energy from celestial sources, not the energy from system noise) from the antennas add together in-phase for maximum gain while the uncorrelated energy (aka system noise) adds together randomly and the correalted energy from elsewhere in the primary beam adds together with sub-optimal phasing that varies across the primary beam. Determining how much gain each primary beam location experiences due to summing is the main objective of the phased array beam pattern analysis.

The time alignemnt depends on the _geometric delays_ from the antennas to the array reference location as the RF wavefront propgates across the array. These geometric delays differ depending on phase center location and array geometry. The signal from each antenna also experiences a finite _fixed delay_ that generally differs between antennas (e.g. due to cable length variations). The geometric delays vary as the phase center moves relative to the array whereas the fixed delays, as their name implies, do not generally vary.  Geometric delays can be calculated when the phase center and array geometry are known, but fixed delays must be measured in situ or given _a priori_.  Other delays, such as those induced by atmosphereic variations over the array, can further complicate things but this analysis will ignore these other sources of delay and delay variation.

Once the geometric and fixed delays are known and compensated for, the signals must still be phase aligned before they will add contructively.  The signal from an antenna can be delay corrected for a single delay value only.  This delay value is the geometric delay for the _delay center_ (generally the same location as the phase center) plus the antenna's fixed delay. Delay corrections are generally performed in discrete steps, which means that there is usually a small residual delay that remains uncorrected.  Phase is the product of frequency and time, ``\phi = \omega \tau``, so the residual delay (a time offset) results in a frequency dependent phase offset.  This phase offset can be removed on a per channel basis after the signal has been channelized into multiple frequency components (e.g. by an FFT or PFB). Locations other than the delay center also have their own residual delays that vary depending on their relationship to the delay center.  Additionally, each antenna has a unique instrumental phase that results from many different factors and must be measured in situ.  Determination of the instrumental phases is beyond the scope of this analysis, so we assume here that the intrumental phases are zero or have already been corrected for.
"""

# ╔═╡ e7853568-21a5-11eb-13de-f17b68612831
md"""Phasing the antennas to a specific phase center requires determination of the geometric delays.  The concept for this is relatively straightforard.  The geometric delay is proportional to the geometric distance that the planar wavefront travels from antenna to array reference location. This distance is the projection of the baseline vector ``\mathbf{b}``, defined to be the vector from the array reference position to the antenna position, onto the unit vector ``\hat{s}``, which points in the direction of the delay center from the array reference position.  This can be computed by taking the dot product of ``\mathbf{b}`` and ``\hat{s}``, which is written as ``\mathbf{b}\cdot\hat{s}``.  The array reference position is usually treated as the origin of the array reference frame so the baseline vector ``\mathbf{b}`` is the same as the antenna position in the array reference frame.  To take the dot product, both ``\mathbf{b}`` and ``\hat{s}`` must be represented in the same frame."""

# ╔═╡ 190282de-21a2-11eb-176e-178fd029cb6d
md"""## Primary Beam Coordinates"""

# ╔═╡ c2295c46-219f-11eb-2d92-cfdeba1d942d
md"""The analysis starts by exploring the coordinates of various points within the primary beam of the relescope. The primary beam is modeled as a circle with an angular diameter equal to the full width half maximum (FWHM) of a single radio telescope. For convenience, this circle is inscribed in a square region of sky which is represented by a two dimensional grid of points with the central point corresponding to the center of the primary beam, i.e. the pointing center. To keep things simple (or at least less complicated), this analysis uses a _flat sky_ approximation.
"""

# ╔═╡ 54a02260-219d-11eb-17b6-17b0e0bc9d8d
md"The axes of this grid correspond to cross-elevation (horizontal, aka `xel`) and elevation (vertical, aka `el`) as seen from the telescope's perspective. The origin of the grid is the center of the primary beam. The unit steps in the cross-elevation and elevations directions are `dxel` and `del`, respectively, though they are sometimes shortened here to `dx` and `de`. The coordinate frame of this grid is referred to here as the _primary beam frame_. A primary beam with FWHM of 2 (arbitrary angular units) and `dxel`, `del` values of 0.25 (arbitrary angular units) can depicted as: 
"

# ╔═╡ bc7cb314-2150-11eb-2900-9b1c7ab0a8d9
begin
	plot(cos, sin, 0:0.01:2π, aspect_ratio=:equal, legend=false, xflip=true)
	scatter!(vec([(x,y) for y=-1:0.25:1, x=-1:0.25:1]), ms=2)
	title!("Primary Beam")
	xlabel!("Cross Elevation (dxel=0.25)")
	ylabel!("Relative Elevation (del=0.25)")
end

# ╔═╡ 2f22727e-22ef-11eb-379f-fd85d1a81d55
md"""For reasons that will be clear later, we have opted to orient the cross elevation (horizontal) axis with poitive values to the left of the origin and negative values to the right."""

# ╔═╡ a018417a-2227-11eb-2f6d-b7460d2e140d
md"""After compensating for geometric delay at the delay center, the residual delay at the delay center is zero.  The rest of the antenna's primary beam will have a delay gradient that depends on the geometry of the delay center relative to the antenna location and array reference location.  To compute the delay for a given point, we must create a unit vector, ``\hat{s}``, in the direction of the point and take the dot product of that and the baseline vector ``\mathbf{b}``.  We are free to choose the coordinate frame in which ``\hat{s}`` is represented.  Usually that frame will be different from the frame in which the antenna positions are expressed, so before taking the dot product we need to ensure we convert the antenna positions into the frame of ``\hat{s}`` or vice versa.  For beam pattern analysis we will typically have more ``\hat{s}`` vectors than antennas, so it will be more efficient to transform the antenna positions into the frame of the ``\hat{s}`` vectors."""

# ╔═╡ 9eb14aca-2323-11eb-3543-33f08ff00e2c
md"""As described above, the coordinate frame for ``\hat{s}`` is most naturally expressed in terms of spherical coordintates of ``(\theta, \phi, 1)``, where ``\theta`` is ``\sphericalangle xel``, ``\phi`` is ``\sphericalangle el``, and ``1`` is the unit length of ``\hat{s}``. The ``\hat{s}`` vector needs to be expressed in a rectalinear form so that we can take the dot product with rectalinear baseline vector ``\mathbf{b}``. The rectalinear form of ``\hat{s}`` is also known as _direction cosines_. The ERFA library we are using provides a spherical to Cartesian function, ``\mathtt{s2c}(θ,ϕ)``, that we can use to calculate ``\hat{s}`` for any values of `xel` and `el`.  Here is how to calculate ``\hat{s}`` for the pointing center and two other points (one on each axis):"""

# ╔═╡ 5b49716c-2327-11eb-2df4-37ceb3707973
begin
	dxel = 0.25 * ERFA.DD2R
	del  = 0.25 * ERFA.DD2R
	[ERFA.s2c(θ, ϕ) for (θ, ϕ) in [(0,0), (dxel,0), (0,del)]]
end

# ╔═╡ a85178f6-2327-11eb-2535-c1f23b6a5484
md"""As you can see, `s2c()` takes two angles (NB: in radians) and returns a rectalinear vector that, in our case, corresponds to `(s, xel, el)`, where `s` is the direction of the source (away from the viewer), `xel` is the cross elevation direction, and `el` is the elevation direction. We choose `xel` to be positive to the left of the origin and `el` positive above the origin to make this a right-handed system because the antenna positions are also specified in a right handed system."""

# ╔═╡ a0fd25d6-2335-11eb-223d-fb1a8be4a076
md"""## Antenna Position Coordinates"""

# ╔═╡ ac2e3c86-222b-11eb-0297-f1aa2b01b60d
md"""Antenna positions are often represented in a topocentric (East, North, Up) or ENU frame, which is a right-handed system.  Antenna positions can also be represented in an _Earth Centered Earth Fixed_ (ECEF) frame known as the _International Terrestrial Reference Frame_ (ITRF).  The ECEF frame can be converted to the ENU frame used in the analysis here, so we will not treat it separately.  Antenna positions are usually given in meters, but other units are sometimes used as well (e.g. nanoseconds or wavelengths at a given frequency).  For our purposes, we will ultimately need to convert the units to time and/or a frequency dependent phase."""

# ╔═╡ 5493615c-2336-11eb-2c36-0772eadafb39
md"""## Coordinate Frame Transformations"""

# ╔═╡ 396e99c8-2336-11eb-24a9-fddf80b160b5
md"""The (s, xel, el) primary beam frame described above is beam-centric in that it does not depend on the location of the pointing center on the sky.  On the other hand, the transformation from the topocentric ENU frame to the beam-centric (s, xel, el) frame does depend on the location of the pointing center."""

# ╔═╡ ff98a0d8-2252-11eb-1056-59c804698c62
md"""The transformation from the ENU frame to the (s, xel, el) frame can be represented as a series of two rotations around various axes of the frame iteslf. Transforming from the ENU frame to the (s, xel, el) frame can be performed by:

  1. Rotating the (E,N,U) frame around the U (Up) axis in the clockwise direction (as viewed looking down to the origin) by the topocentric azimuth angle of the pointing center minus 90 degrees.  This results in an intermediate (E',N',U) frame.
  2. Rotating the intermediate (E',N',U) frame around the N' axis in the clockwise direction (as viewed looking at the origin from positive N') by the elevation angle.  This results in a (E", N', U') frame that is equivalent to the (s, xel, el) frame in which our ``\hat{s}`` vectors are expressed.
"""

# ╔═╡ f4f6ab3a-222c-11eb-0564-17ceb3dde858
md"""Each of these rotations can be represented by a 3x3 rotation matrix. These matices can be multiplied together, in the proper order, to produce a single composite 3x3 rotation matrix that transforms the ENU frame into the beam-centric (s, xel, el) frame for the azimuth and elevation of a given the pointing center.  The ENU antenna position(s) represented as a 3 element column vector (or 3xM matrix of M antenna positions) can be left multiplied by this composite rotation matrix to produce the antenna position(s) in the (s, xel, el) frame.  Here we define an `enu2beam(az, el)` function that returns this composite rotation matrix for a pointing center at azimuth `az` and elevation `el` (both given in degrees):"""

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
md"""Notice that in the last example (zenith from south), the sign of `xel` and `el` will be opposite from the example that preceeds it (zenith from north) even though both have the local zenith as the pointing center."""

# ╔═╡ e8db8432-2333-11eb-3195-877803e736ba
md"""## Putting It All Together"""

# ╔═╡ e8b55f50-2333-11eb-0f6c-abc07cdfd3cc
md"""We now have the tools to start building and putting together all the various parts of the analysis:
  - Grid of points covering the primary beam
  - Unit ``\hat{s}`` vectors for each of those points
  - Antenna positions in ENU form
  - Conversion from ENU frame to beam frame for arbitrary pointing center
  - Dot product to compute distance/delay from antenna to array reference poistion
"""

# ╔═╡ d620f9aa-23ba-11eb-11a2-75d99569a381
md"""Here we create a grid of points to cover the primary beam and unit ``\hat{s}`` vectors corresponding to those points:"""

# ╔═╡ 165f5502-23b6-11eb-3824-dde609547349
npoints = 65

# ╔═╡ 5a17d3a8-23b6-11eb-278f-0f399b01b6b1
fwhm = 2.0

# ╔═╡ 11735a72-23b7-11eb-0c75-39f94e10ba04
xels = range(-fwhm/2, +fwhm/2, length=npoints)

# ╔═╡ 4697cad0-23b7-11eb-17fd-d73044fab8bb
els = range(-fwhm/2, +fwhm/2, length=npoints)

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
md"""With the antenna position in the same frame as the unit ``\hat{s}`` vectors, we can take the dot products to calculate the distances from the antenna to the array reference position (i.e. the origin of the ENU frame) for each point in the primary beam.  We are more interested in the relative distances across the beam rather than the absoulte distances.  To calculate the relative distances we subtract each unit vector ``\hat{s}`` from unit vector ``\hat{d}`` that points in the direction of the delay center (assumed here to be the same as the pointing center).  These relative distances can be plotted as a contour plot.
"""

# ╔═╡ a2bf4b28-23d5-11eb-1ef1-1170b39fcb9e
reldists = map(s->dot(ERFA.s2c(0,0)-s, user_antpos_pointing), grid_svecs);

# ╔═╡ 601ffcd6-23c4-11eb-2012-01bb72855b38
contour(xels, els, reldists',
	aspect_ratio=:equal, fill=true, clabels=true,
	fillcolor=:lightrainbow, linecolor=:black,
	title="Antenna to Origin: Distance\nrelative to (az=$(user_az), el=$(user_el))",
	xlabel="Cross Elevation (degrees)",
	ylabel="Elevation (degrees)"
)

# ╔═╡ a37ac5d2-2483-11eb-3761-c18b249d25c1
md"""For a given frequency, these relative distances can be converted to relative phases which can also be plotted as a heapmap plot:

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
	plt=heatmap(xels, els, relphases',
		aspect_ratio=:equal, color=:hsv,
		clims=(-180,180), levels=71,
		title="Antenna to Origin: Degrees of Phase\nrelative to az=$(user_az), el=$(user_el) at $(user_mhz) MHz (λ=$(round(1e5λ)/1e3) cm)\n",
		xlabel="Cross Elevation (degrees)",
		ylabel="Elevation (degrees)"
	)
	# Plot 5 degree contours if pases are within +/- 90 degrees
	# (contouring around phase wraps gets ugly)
	if maximum(abs.(relphases)) <= 90
		contour!(xels, els, relphases',
			aspect_ratio=:equal, fill=false, clabels=true,
			levels=71, linecolor=:black,
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
md"""## Adding a Second Antenna"""

# ╔═╡ 4196dd3a-24b3-11eb-2e0d-650a39dc4bb6
md"""Now that we can compute phases across the primary beam of an antenna relative to the array center for a given pointing/delay center at a given frequency, the next thing to do is add another antenna.  With two antennas, we can delay correct them for a common delay/pointing center, calcluate their phases across the primary beam, convert these phases to unit vectors in the complex plane (i.e. complex numbers), and add corresponding complex numbers together to create the complex beam pattern of the coeherent beam created by adding the two delay corrected signals together."""

# ╔═╡ aef63f92-24b3-11eb-1a6c-79b47130a4a9
md"""To keep the UI simple, we will use new antennas, new pointing center, and new frequency.  The new antennas will have an `Up` coordinate of zero, but you can still move them East/West and North/South."""

# ╔═╡ dc9d05c2-24b8-11eb-0ec2-87d8c461845a
md"""
$(@bind ant1_east html"<input type=range min=-100000 max=100000 value=3000 style='width: 75%'> Ant 1 East</br>")
$(@bind ant1_north html"<input type=range min=-100000 max=100000 value=0 style='width: 75%'> Ant 1 North</br>")
$(@bind ant2_east html"<input type=range min=-100000 max=100000 value=0 style='width: 75%'> Ant 2 East</br>")
$(@bind ant2_north html"<input type=range min=-100000 max=100000 value=3000 style='width: 75%'> Ant 2 North</br>")
$(@bind ant12_az html"<input type=range min=0 max=359 value=90 style='width: 75%'> Ants 1,2 Azimuth</br>")
$(@bind ant12_el html"<input type=range min=0 max=90 value=0 style='width: 75%'> Ants 1,2 Elevation")
$(@bind ant12_mhz html"<input type=range min=500 max=2000 value=1420 style='width: 75%'> Ants 1,2 Frequency")
"""

# ╔═╡ 15fea64e-24ba-11eb-09cf-4d47980c7fb5
begin
	ant12_enu2beam = enu2beam(ant12_az, ant12_el)
	
	ant1_antpos_pointing = ant12_enu2beam * [ant1_east, ant1_north, 0] ./ 100
	ant2_antpos_pointing = ant12_enu2beam * [ant2_east, ant2_north, 0] ./ 100
	
	ant1_reldists = map(s->dot(ERFA.s2c(0,0)-s, ant1_antpos_pointing), grid_svecs)
	ant2_reldists = map(s->dot(ERFA.s2c(0,0)-s, ant2_antpos_pointing), grid_svecs)
	
	λ12 = ERFA.CMPS / ant12_mhz / 1e6
	ant1_relphases = rem2pi.(2π .* ant1_reldists ./ λ12, RoundNearest)
	ant2_relphases = rem2pi.(2π .* ant2_reldists ./ λ12, RoundNearest)
	#ant1_relphases = 2π .* ant1_reldists ./ λ12
	#ant2_relphases = 2π .* ant2_reldists ./ λ12

	
	ant1z = exp.(1im*ant1_relphases)
	ant2z = exp.(1im*ant2_relphases)
	ant12z = ant1z .+ ant2z
	
	p1=heatmap(xels, els, ant1_relphases' .* ERFA.DR2D,
		aspect_ratio=:equal, color=:hsv,
		clims=(-180,180), levels=71,
		title="Ant 1 Phases",
		#xlabel="Cross Elevation (degrees)",
		#ylabel="Elevation (degrees)"
	)

	p2=heatmap(xels, els, ant2_relphases' .* ERFA.DR2D,
		aspect_ratio=:equal, color=:hsv,
		clims=(-180,180), levels=71,
		title="Ant 2 Phases",
		#xlabel="Cross Elevation (degrees)",
		#ylabel="Elevation (degrees)"
	)
	
	p3=heatmap(xels, els, abs.(ant12z)',
		aspect_ratio=:equal, fill=true, fillcolor=:thermal,
		climes=(0,2), levels=9,
		title="Coherent Beam (amplitude)"
	)

	p4=heatmap(xels, els, angle.(ant12z)' .* ERFA.DR2D,
		aspect_ratio=:equal, color=:hsv,
		clims=(-180,180), levels=71,
		title="Coherent Beam (phase)"
	)

	plot(p1, p2, p3, p4)
end

# ╔═╡ 63963368-24bf-11eb-1237-773f333dcb96
ant1_enu = [ant1_east, ant1_north, 0] ./ 100

# ╔═╡ d853dd02-24bf-11eb-11f2-6d197e32c0c3
ant2_enu = [ant2_east, ant2_north, 0] ./ 100

# ╔═╡ e41097de-24bf-11eb-2fed-a527d5c53dbe
ant12_azel = (ant12_az, ant12_el)

# ╔═╡ 07431902-24c0-11eb-2f8c-07c42874be62
ant12_freq_mhz = ant12_mhz

# ╔═╡ b9e16b1a-23ba-11eb-08ad-8fc21f98dca9
md"""---"""

# ╔═╡ 3beff902-225b-11eb-267d-21ef072bf855
ant1_reldists

# ╔═╡ 428408b8-2151-11eb-0530-d1996a93b332


# ╔═╡ 3e4354d4-2151-11eb-384f-7706d85ce220


# ╔═╡ bccfe344-214c-11eb-2a1a-ff912a978a1e


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
# ╟─5a17d3a8-23b6-11eb-278f-0f399b01b6b1
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
# ╠═de5be0d6-2489-11eb-21b1-e3218333e30d
# ╟─7324d316-24ae-11eb-2112-a51606afab74
# ╟─22338236-24b3-11eb-02aa-0bc76375dfd3
# ╟─4196dd3a-24b3-11eb-2e0d-650a39dc4bb6
# ╟─aef63f92-24b3-11eb-1a6c-79b47130a4a9
# ╟─dc9d05c2-24b8-11eb-0ec2-87d8c461845a
# ╟─15fea64e-24ba-11eb-09cf-4d47980c7fb5
# ╟─63963368-24bf-11eb-1237-773f333dcb96
# ╟─d853dd02-24bf-11eb-11f2-6d197e32c0c3
# ╟─e41097de-24bf-11eb-2fed-a527d5c53dbe
# ╟─07431902-24c0-11eb-2f8c-07c42874be62
# ╠═b9e16b1a-23ba-11eb-08ad-8fc21f98dca9
# ╠═3beff902-225b-11eb-267d-21ef072bf855
# ╠═428408b8-2151-11eb-0530-d1996a93b332
# ╠═3e4354d4-2151-11eb-384f-7706d85ce220
# ╠═bccfe344-214c-11eb-2a1a-ff912a978a1e
