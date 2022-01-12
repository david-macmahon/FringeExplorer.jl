### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ 89001295-1da7-496b-8c72-475520512b8a
begin
	import Pkg
	Pkg.activate(Base.current_project())
	using CairoMakie
	using GeometryBasics
	using PlutoUI
	using Rotations
	using FringeExplorer
	using LinearAlgebra
	using ERFA # For sepp until we add a composable version to FringeExplorer.
	md"""# Coordinate Explorations"""
end

# ╔═╡ 3639b720-f3c6-49e9-9106-84398cc332df
begin
	f=Figure()
	a=Axis(f[1,1], limits=((-5,5),(-5,5)), aspect=1)
	dy = sin(π/3)
	for y in -4:4
		xmin = -4+abs(y/2)
		xmax = - xmin
		for x in xmin:xmax
			lines!(Circle(Point2f0(x,y*dy), 0.5))
		end
	end
	f
end

# ╔═╡ e1e798aa-259a-435b-baea-857bf893cde3
begin
	f2=Figure()
	a2=Axis(f2[1,1], limits=((-5,5),(-5,5)), aspect=1)
	lines!(Circle(Point2f(0), 0.5))
	for r in 1:4
		nx = 6r
		for i in 0:nx-1
			θ = i*2π/nx
			x = r*cos(θ)
			y = r*sin(θ)
			lines!(Circle(Point2f(x,y), 0.5))
		end
	end
	f2
end

# ╔═╡ 3319a555-62b0-44ac-a2e9-01ebeca5a07b
fwhm = deg2rad(4)

# ╔═╡ 62fa9f61-3aa4-4542-9629-67e93198d5b1
α0=5π/4

# ╔═╡ 48266d0e-c873-4597-a981-e86c47250a8b
δ0=π/6

# ╔═╡ 34b142aa-0b8f-4311-b922-88fddc7fdc15
centers=beamrings(α0, δ0; nrings=4, dϕ=fwhm)

# ╔═╡ f657d533-e777-4ec0-8367-c035279d9d20
centerc = FringeExplorer.s2c(centers)

# ╔═╡ 10fb9c7d-bcc3-45ce-8cf2-eca33ab3a5d0
# Cross beam center vectors with Z axis to get rotation axes
crosses = cross.(Ref([0.0,0.0,1.0]), eachcol(centerc))

# ╔═╡ 0417c07f-47d3-4c12-8960-333c78f3288a
beamseps = sepp.(Ref([0.0,0.0,1.0]), eachcol(centerc))

# ╔═╡ 8c965934-14da-4bf8-b30c-8ad94ee601cd
begin
	f3=Figure()
	a3=Axis3(f3[1,1], limits=((-1,1),(-1,1),(-1,1)),
		              aspect=(1,1,1),
					  azimuth=deg2rad(210))
	s3=Sphere(Point3f(0), 1.0)
	c3=Circle(Point2f(0), 1.0)#tan(fwhm/2))
	#arrows!([0,0,0], [0,0,0], [0,0,0], [1,0,0], [0,1,0], [0,0,1])
	mesh!(s3, color=RGBAf(1.0, 1.0, 1.0, 0.4))
	lines!(c3, color=RGBAf(0.5, 0.5, 0.5, 0.5))
	for ax in (Vec3f(0,1,0), Vec3f(1,0,0))
		transformation=Transformation(;rotation=qrotation(ax, π/2))
		lines!(c3; color=RGBAf(0.5, 0.5, 0.5, 0.5), transformation)
	end
	for (i, beam) in enumerate(eachcol(centerc))
		trans = Transformation(;scale=Vec3f(tan(fwhm/2)),
								translation=beam,
								rotation=qrotation(Vec3f(crosses[i]), beamseps[i]))
		lines!(c3, transformation=trans)
	end
	f3
end
	

# ╔═╡ Cell order:
# ╟─89001295-1da7-496b-8c72-475520512b8a
# ╟─3639b720-f3c6-49e9-9106-84398cc332df
# ╟─e1e798aa-259a-435b-baea-857bf893cde3
# ╠═3319a555-62b0-44ac-a2e9-01ebeca5a07b
# ╠═62fa9f61-3aa4-4542-9629-67e93198d5b1
# ╠═48266d0e-c873-4597-a981-e86c47250a8b
# ╠═34b142aa-0b8f-4311-b922-88fddc7fdc15
# ╠═f657d533-e777-4ec0-8367-c035279d9d20
# ╠═10fb9c7d-bcc3-45ce-8cf2-eca33ab3a5d0
# ╠═0417c07f-47d3-4c12-8960-333c78f3288a
# ╠═8c965934-14da-4bf8-b30c-8ad94ee601cd
