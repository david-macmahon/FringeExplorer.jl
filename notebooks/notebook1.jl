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

# ╔═╡ 6c692a11-cb4c-425d-9acb-9529bd382dd2
begin
    import Pkg
	Pkg.activate(Base.current_project())

	using FringeExplorer, Latexify, RadioInterferometry, Symbolics, Plots, PlotlyJS, LaTeXStrings, PlutoUI
	Latexify.set_default(fmt=FancyNumberFormatter(4, "{\\cdot}"))
	#Latexify.set_default(fmt=FancyNumberFormatter("%.1g",
	#	SubstitutionString{String}("\\g<mantissa> {\\\\cdot} 10^{\\g<exp>}")))
	#LaTeXString(numformat(latexify(FringeExplorer.wdw_m)))
	md"""
	# Fringe Explorer
	This notebook is for exploring the relative geometric delays, and their rates of change, as seen by antennas in a terrestrial radio telescope array as it tracks a sidereal source moving across the sky.

	Of course we now know that the ground based antennas are doing (most of) the moving thanks to Earth's rotation, but as observers who live in the same reference frame as the antennas it can sometimes be convenient to think in terms of the sky moving instead.
	"""
end

# ╔═╡ 655091ad-69a0-4672-abb4-37ee03b27cff
md"""
## Getting Oriented
### XYZ frame
For the purposes of this discussion we will represent antenna positions, and the array's reference position, in a right-handed XYZ frame where the Z axis is parallel to the rotational axis of Earth and the XY plane parallel to the equarorial plane of Earth.  One such frame commonly used for antenna positions is the *International Terrestrial Reference Frame* (ITRF), a geocentric frame where Z is the rotational axis of Earth, X is in the equatorial lane and passes through the prime meridian (i.e. longitude 0), and Y completes the right hand frame (i.e. passing though longitude 90 degrees).  The frame used here is a modified version of this frame where the origin has been translated to the array reference position and the frame has been rotated around the Z axis such that X points to the longitude of the array reference position.
"""

# ╔═╡ 7ff31ce8-190e-4018-ba49-c2f11ed82801
md"""
### Hour Angle and Declination frame
We will also choose to represent sidereal positions as an *hour angle*, represented by the variable ``h``, and *declination* angle, represented by the variable ``\delta``.  Being angles, we will store and calculate with them in radians, but often display them as degrees.  This choice of representation is convenient because for a given sidereal position, ``h`` is a linear function of time and ``\delta`` is a constant.  Positions on the local meridian have ``h = 0``, positions to the west have ``h > 0``.  Using ``\Omega`` to represent rate of rotation, ``h`` for a given sideral source can be expressed as ``h = \Omega t``, with ``t = 0`` corresponding to the time when the source is *transiting*, i.e. on/crossing the local meridian.
"""

# ╔═╡ 423f2ef7-932d-4eb8-8ea6-304e2b7cbf2e
md"""
### UVW frame
With these conventions established, we can now discuss how to compute the relative geometric delays between wavefront arrival times at a given antenna and the array reference position.  One way to do this is to rotate the XYZ frame we are using for antenna positions into a new frame known as the UVW frame (also right handed).  The U axis of the UVW frame points eastward, the V axis points northward, and the W axis points at the source.  Because the W axis points to the source that is moving, the UVW frame is not a fixed frame, but a rotating frame.  When representing the antenna positions in the UVW frame, the W component of this position vector represents to *line of sight* distance between the antenna and the reference position as seen from the source position.  In other words, the RF wavefront arrives at one end of the vector first, then travels a distance W before arriving at the other end of the vector.  For our purposes here we will ignore the amount of Earth rotation that happens during this time.  The XYZ coordinates of the antennas are typically specified in meters, but speed of light can be used to convert W from meters to seconds, or more typically nanoseconds, of delay which can be used along with observing frequency to convert W to units of wavelength.
"""

# ╔═╡ c02e8687-be77-4af0-a0a1-6901ad02f85e
md"""
Rotating antenna position vectors from the XYZ frame to the UVW frame can be accomplished by rotating our XYZ frame about the Z (i.e. third) axis by ``h``, producing a (X',U,Z) frame, then rotating that frame about the U (i.e. second) axis by ``\delta``, producing a (W,U,V) frame which can then permuted to (U,V,W) where U is east, V is north, and W is in the direction of ``h`` and ``\delta``.  This operation can be performed by left multiplying the antenna position vector by a 3x3 rotation matrix.

The `RadioInterferometry.jl` package provides a `xyz2uvw` function that can generate these 3x3 rotation matrices for any given ``h`` and ``\delta``.

!!! tip
    The `xyz2uvw` function supports XYZ frames that, unlike ours, are oriented with X parallel to the great circle containing the prime meridian (e.g. the ITRF XYZ frame).  This involves adjusting ``h`` by the logitude of the array's reference position, but since our X axis is already oriented to the local meridian we can simply pass 0 for longitude.

The rotation from XYZ to UVW for a given hour angle ``h`` and declination ``\delta`` can be expressed as:
"""

# ╔═╡ 0c19d2b3-3869-48f1-b2eb-1ff157e931cf
begin
	@variables x y z u v w h δ
	([u,v,w] ~ :($(xyz2uvw(h, δ, 0)) * $([x,y,z])))|>latexify
end

# ╔═╡ 11567a37-55b2-4d9e-88fb-9efdf7ada485
md"""
## Play time!

Change the sliders here to select between hour angle (in minutes) and declination (in degrees).  The values you selected are shown below followed by the XYZ to UVW transformation for the given parameters.

$(@bind h0 Slider(-360:5:360, default=0)) Hour angle (minutes)

$(@bind δ0 Slider(-90:5:90, default=0)) Declination (degrees)
"""

# ╔═╡ 1ceb4ee8-2bbd-45a7-b5a4-886c949c17e3
(h=h0, δ=δ0)

# ╔═╡ aee119f3-b190-4342-9d48-c0ee9f8e2859
begin
	([u,v,w] ~ :($(xyz2uvw(deg2rad(h0/4), deg2rad(δ0), 0)) * $([x,y,z])))|>latexify
end

# ╔═╡ 979c1962-db9a-4581-95a8-15b183afd990


# ╔═╡ d4145813-5e18-4230-8389-d8935d72ba4b
FringeExplorer.wdw_m

# ╔═╡ Cell order:
# ╟─6c692a11-cb4c-425d-9acb-9529bd382dd2
# ╟─655091ad-69a0-4672-abb4-37ee03b27cff
# ╟─7ff31ce8-190e-4018-ba49-c2f11ed82801
# ╟─423f2ef7-932d-4eb8-8ea6-304e2b7cbf2e
# ╟─c02e8687-be77-4af0-a0a1-6901ad02f85e
# ╟─0c19d2b3-3869-48f1-b2eb-1ff157e931cf
# ╟─11567a37-55b2-4d9e-88fb-9efdf7ada485
# ╟─1ceb4ee8-2bbd-45a7-b5a4-886c949c17e3
# ╟─aee119f3-b190-4342-9d48-c0ee9f8e2859
# ╠═979c1962-db9a-4581-95a8-15b183afd990
# ╠═d4145813-5e18-4230-8389-d8935d72ba4b
