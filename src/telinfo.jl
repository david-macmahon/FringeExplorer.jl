# Functions related to telinfo.yml files

using YAML: load_file

"""
    telinfo(yamlfile)

Loads YAML file `yamlfile` using `dicttype=Dict{Symbol,Amy}`.
"""
telinfo(fname) = YAML.load_file(fname, dicttype=Dict{Symbol,Any})

"""
    antpos(telinfo_yml::AbstracyString)
    antpos(telinfo::AbstractDict{Symbol, Any})
    antpos(antennas::AbstractVector{<:AbstractDict{Symbol,Amy}})

Return antenna positions from `antennas` Vector as a 3xN matrix (i.e. one
antenna per columm).  The antenna ordering is sorted by the `number` of each
antenna.

# Examples

```julia
ti = telinfo("telinfo.yml") # => Dict{Symbol, Any}

ap = antpos("telinfo.yml")  # => Array{Float64,2}
ap = antpos(ti)             # => Array{Float64,2}
ap = antpos(ti[:antennas])  # => Array{Float64,2}
```

!!! note
    The antenna numbers are used for sorting, but they are NOT antenna indexes.
"""
function antpos(antennas::AbstractVector{<:AbstractDict{Symbol,Any}})
    ants = sort(antennas, by=a->a[:number])
    reduce(hcat, (a[:position] for a in ants))
end,

function antpos(ti::AbstractDict{Symbol,Any})
    antpos(ti[:antennas])
end,

function antpos(tiyml::AbstractString)
    antpos(telinfo(tiyml))
end
