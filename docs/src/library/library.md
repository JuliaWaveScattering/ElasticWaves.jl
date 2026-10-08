# Library

```@meta
CurrentModule = ElasticWaves
```

All the types and functions which have documentation are listed below.

```@contents
Pages = ["library.md"]
Depth = 2:2
```

## Elastic media and waves

```@autodocs
Modules = [ElasticWaves]
Order   = [:type, :function]
Pages   = ["src/types.jl", "src/elasticity.jl", "src/fields.jl"]
```

## Roller bearings

Types and functions to describe a roller bearing and the data on its boundaries.

```@autodocs
Modules = [ElasticWaves]
Order   = [:type, :function]
Pages   = ["src/bearing.jl", "src/loading-profile.jl", "src/signal_processing.jl"]
```

## Waves in a raceway

Functions to solve for, and to evaluate, the waves in the raceway of a roller bearing.

```@autodocs
Modules = [ElasticWaves]
Order   = [:type, :function]
Pages   = ["src/cylindrical/elastic_wave.jl", "src/cylindrical/fields.jl", "src/cylindrical/t-matrix.jl"]
```

## Scattering from spheres

```@autodocs
Modules = [ElasticWaves]
Order   = [:type, :function]
Pages   = ["src/source.jl", "src/spherical/fields.jl", "src/spherical/t-matrix.jl"]
```
