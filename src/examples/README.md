# Examples and Driver Scripts

This directory contains example scripts and driver programs for the iFEM-Maxwell-PML project.

## Phase 1: Geometry Generation

### Quick Start

```matlab
% Run the main driver
run('src/examples/cube_sphAir_sphPML_driver.m')

% Or run the test script
run('src/examples/test_geometry.m')
```

### Files

- `cube_sphAir_sphPML_driver.m` - Main driver for Phase 1 geometry testing
- `test_geometry.m` - Comprehensive test script for geometry validation

### Usage

The main driver script demonstrates:

1. **Geometry Generation**: Creates three-region mesh (SiC cube + air sphere + PML shell)
2. **Region Tagging**: Hybrid approach with boundary detection
3. **PML Grading**: Smooth mesh size variation in PML region
4. **Boundary Detection**: Identifies SBC boundary faces
5. **Validation**: Checks geometry correctness and mesh quality

### Parameters

- **SiC Cube**: L = 1 μm (side length)
- **Air Sphere**: r = 6 μm (radius)
- **PML Shell**: r = 8 μm (outer radius, thickness = 2 μm)
- **Wavelength**: λ = 5 μm
- **Mesh Sizes**: h_air = 0.5 μm, h_SiC = 0.25 μm, h_out = 1.0 μm

### Output

- `node` - Node coordinates [N × 3]
- `elem` - Element connectivity [M × 4] (tetrahedra)
- `regionID` - Region ID for each element [M × 1]
- `bdFlag` - Boundary face flags [F × 1]
- `meta` - Metadata structure with counts and statistics

### Visualization

The driver includes automatic visualization showing:
- Region coloring (SiC = red, Air = blue, PML = green)
- Mesh structure and quality
- SBC boundary faces

### Dependencies

- MATLAB's `distmesh3d` (if available)
- `uniformrefine3` from iFEM
- Custom geometry functions in `src/geometry/`

### Troubleshooting

If `distmesh3d` fails, the code will fall back to a simple structured mesh. For production use, consider using GMSH or other professional meshers.

### Next Steps

After successful geometry generation, proceed to:
- Phase 2: Scattering Boundary Condition (SBC)
- Phase 3: Perfectly Matched Layer (PML)
- Phase 4: Post-processing
