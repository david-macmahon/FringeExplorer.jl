# Functions related to telinfo.yml files

import YAML: load_file
using TOML: parse

"""
    telinfo(telinfo_file; dicttype=Dict{Symbol,Any})

Loads telescope information from `telinfo_file`.  If `telinfo_file` ends with
`.yml` or `.yaml` it will be loaded as a YAML file using `dicttype`
dictionaries.  If `telinfo_file` ends with `.toml` it will be loaded as a TOML
file and `dicttype` will be ignored since TOML files always use
`Dict{String,Any}`.
"""
function telinfo(fname; dicttype=Dict{Symbol,Any})
    _, ext = splitext(fname)
    if ext == ".yml" || ext == ".yaml"
        load_file(fname; dicttype)
    elseif ext == ".toml"
        parse(fname)
    else
        error("""unknown extention: "$(ext)\"""")
    end
end

"""
    antpos(telinfo_yml::AbstracyString)
    antpos(telinfo::AbstractDict{Symbol, Any})
    antpos(telinfo::AbstractDict{String, Any})
    antpos(antennas::AbstractVector{<:AbstractDict{Symbol,Any}})
    antpos(antennas::AbstractVector{<:AbstractDict{String,Any}})

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
end

function antpos(antennas::AbstractVector{<:AbstractDict{String,Any}})
    ants = sort(antennas, by=a->a["number"])
    reduce(hcat, (a["position"] for a in ants))
end

function antpos(ti::AbstractDict{Symbol,Any})
    antpos(ti[:antennas])
end

function antpos(ti::AbstractDict{String,Any})
    antpos(ti["antennas"])
end

function antpos(tiyml::AbstractString; dicttype=Dict{Symbol,Any})
    antpos(telinfo(tiyml; dicttype))
end
