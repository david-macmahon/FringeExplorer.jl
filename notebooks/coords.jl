### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ 89001295-1da7-496b-8c72-475520512b8a
begin
	import Pkg
	Pkg.activate(Base.current_project())
	#using WGLMakie
	using CairoMakie
	using GeometryBasics
	using PlutoUI
	using Rotations
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

# ╔═╡ 25b138a9-e15b-4c52-9a2a-9292d205ad07
begin
	function c2s(x,y,z)
		θ=atan(y,x)
		ϕ=atan(z, hypot(x,y))
		[θ, ϕ]
	end
	function c2s(xyz::AbstractVector{<:Real})
		c2s(xyz[1], xyz[2], xyz[3])
	end
	function c2s(xyz::AbstractMatrix{<:Real})
		#reduce(hcat, c2s.(xyz[1,:], xyz[2,:], xyz[3,:]))
		mapreduce(c2s, hcat, eachcol(xyz))
	end
	(0,1,0)=>c2s(0,1,0)
end

# ╔═╡ a8da9102-89a1-40b6-aebc-64e81286cf4d
begin
	function s2c(θ,ϕ)
		z, r = sincos(ϕ)
		y, x = sincos(θ) .* r
		[x, y, z]
	end
	#function s2c(θϕ::AbstractVector{<:Real})
	#	s2c(θϕ[1], θϕ[2])
	#end
	function s2c(θϕ::AbstractMatrix{<:Real})
		#mapreduce(s2c, hcat, eachcol(θϕ))
		@views reduce(hcat, s2c.(θϕ[1,:], θϕ[2,:]))
	end
	(π/2, 0) => s2c(π/2, 0)
end

# ╔═╡ 89058f56-f2db-4aa9-8ce7-f0433e889320
function beamrings(h, δ, nrings=4, dϕ=deg2rad(10/3600))
	centers = [[0, π/2]]
	for r in 1:nrings
		nθ = 6r
		for i in 0:nθ-1
			θ = i*2π/nθ
			ϕ = π/2 - r*dϕ
			push!(centers, [θ, ϕ])
		end
	end
	mapreduce(c2s, hcat, Ref(RotZY(h, π/2-δ)) .* s2c.(centers))
end	

# ╔═╡ 19c5dccb-5304-4a40-8c51-3e824e4fe90a
function beamrings2(h, δ, nrings=4, dϕ=deg2rad(10/3600))
	centers = [[0, π/2]]
	for r in 1:nrings
		nθ = 6r
		for i in 0:nθ-1
			θ = i*2π/nθ
			ϕ = π/2 - r*dϕ
			push!(centers, [θ, ϕ])
		end
	end
	RotZY(h, π/2-δ) * mapreduce(s2c, hcat, centers) |> c2s
end	

# ╔═╡ 6142df3c-70f6-47c9-96e1-26ac1e7afe60
function beamrings3(θ, ϕ, nrings=4, dϕ=deg2rad(10/3600))
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

# ╔═╡ eceac8ca-ce9b-42b1-8569-fc6b4cccf324
beamrings(0, 1)

# ╔═╡ 79c852fa-ce49-4b31-a615-086ac927c4d2
beamrings2(0, 1)

# ╔═╡ 34b142aa-0b8f-4311-b922-88fddc7fdc15
beamrings3(0,1)

# ╔═╡ Cell order:
# ╠═89001295-1da7-496b-8c72-475520512b8a
# ╟─3639b720-f3c6-49e9-9106-84398cc332df
# ╠═e1e798aa-259a-435b-baea-857bf893cde3
# ╟─25b138a9-e15b-4c52-9a2a-9292d205ad07
# ╠═a8da9102-89a1-40b6-aebc-64e81286cf4d
# ╟─89058f56-f2db-4aa9-8ce7-f0433e889320
# ╟─19c5dccb-5304-4a40-8c51-3e824e4fe90a
# ╠═6142df3c-70f6-47c9-96e1-26ac1e7afe60
# ╠═eceac8ca-ce9b-42b1-8569-fc6b4cccf324
# ╠═79c852fa-ce49-4b31-a615-086ac927c4d2
# ╠═34b142aa-0b8f-4311-b922-88fddc7fdc15
