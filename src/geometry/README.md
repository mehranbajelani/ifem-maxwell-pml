# Geometry Module

This module contains functions for creating and managing the three-region geometry required for electromagnetic scattering simulations.

## Files

- `makeSphericalAirAndPML.m` - Main function for creating the three-region mesh

## Geometry Description

The geometry consists of three nested regions:

1. **SiC Cube** (Region ID = 1)
   - Centered at origin
   - Side length: L = 1 μm
   - Material: Silicon Carbide (SiC)

2. **Air Sphere** (Region ID = 2)
   - Radius: r_SBC = 6 μm
   - Material: Air (ε_r = 1, μ_r = 1)
   - Contains the SiC cube

3. **PML Shell** (Region ID = 3)
   - Inner radius: r_SBC = 6 μm
   - Outer radius: r_PML = 8 μm
   - Thickness: t_PML = 2 μm
   - Material: Perfectly Matched Layer

## Usage

```matlab
% Basic usage
[node, elem, regionID, boundaryID] = makeSphericalAirAndPML();

% With custom parameters
[node, elem, regionID, boundaryID] = makeSphericalAirAndPML(...
    'L', 1e-6, 'r_SBC', 6e-6, 'r_PML', 8e-6, ...
    'h_air', 0.5e-6, 'visualize', true);
```

## Parameters

- `L` - SiC cube side length (default: 1e-6 m)
- `r_SBC` - Air sphere radius (default: 6e-6 m)
- `r_PML` - PML outer radius (default: 8e-6 m)
- `h_air` - Air region mesh size (default: 0.5e-6 m)
- `h_SiC` - SiC region mesh size (default: 0.25e-6 m)
- `h_PML` - PML region mesh size (default: 0.5e-6 m)
- `refine` - Number of mesh refinements (default: 0)
- `visualize` - Show mesh plot (default: false)

## Output

- `node` - Node coordinates [N x 3]
- `elem` - Element connectivity [M x 4] (tetrahedra)
- `regionID` - Region ID for each element [M x 1]
- `boundaryID` - Boundary ID for each face [F x 1]

## Mesh Quality Guidelines

- **Air region**: h_air ≈ λ/10 = 0.5 μm (for λ = 5 μm)
- **SiC region**: h_SiC ≈ 0.25 μm (finer mesh for accuracy)
- **PML region**: h_PML ≈ 0.5 μm (can be coarser)

## Dependencies

- MATLAB's `distmesh3d` (if available)
- `uniformrefine3` from iFEM
- Custom fallback mesh generation

## Testing

Run `test_geometry.m` to validate the geometry generation and mesh quality.
