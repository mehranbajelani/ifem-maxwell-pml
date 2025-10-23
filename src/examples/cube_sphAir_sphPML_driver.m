%CUBE_SPHAIR_SPHPML_DRIVER Main driver for Phase 1 geometry testing
%
% This script demonstrates the usage of makeSphericalAirAndPML.m
% and validates the geometry generation for the electromagnetic
% scattering problem.
%
% Author: iFEM-Maxwell-PML Project
% Date: 2025

clear; clc; close all;

% Setup paths
fprintf('Setting up paths...\n');
% Get the project root (go up one level from src/examples)
project_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(project_root);
addpath(fullfile(project_root, 'src'));
addpath(fullfile(project_root, 'src', 'geometry'));
addpath(fullfile(project_root, 'src', 'physics'));
addpath(fullfile(project_root, 'src', 'post'));
addpath(fullfile(project_root, 'src', 'examples'));

fprintf('=== iFEM-Maxwell-PML: Phase 1 Geometry Driver ===\n\n');

% Physical parameters (SI units)
L = 1e-6;        % SiC cube side length: 1 μm
r_SBC = 6e-6;    % Air sphere radius: 6 μm  
r_PML = 8e-6;    % PML outer radius: 8 μm
lambda = 5e-6;   % Wavelength: 5 μm

% Mesh parameters
h_air = lambda/10;    % Air mesh size: λ/10 = 0.5 μm
h_SiC = L/4;          % SiC mesh size: L/4 = 0.25 μm
h_out = 1.0e-6;       % PML outer mesh size: 1.0 μm

% PML grading parameters
alpha = 1.0;           % PML grading parameter
beta = 1.0;           % PML grading parameter

fprintf('Physical Parameters:\n');
fprintf('  SiC cube: L = %.1f μm\n', L*1e6);
fprintf('  Air sphere: r = %.1f μm\n', r_SBC*1e6);
fprintf('  PML shell: r = %.1f μm (thickness = %.1f μm)\n', r_PML*1e6, (r_PML-r_SBC)*1e6);
fprintf('  Wavelength: λ = %.1f μm\n', lambda*1e6);
fprintf('  Frequency: f = %.1f THz\n', 3e8/lambda/1e12);

fprintf('\nMesh Parameters:\n');
fprintf('  Air mesh size: h_air = %.2f μm (λ/10)\n', h_air*1e6);
fprintf('  SiC mesh size: h_SiC = %.2f μm (L/4)\n', h_SiC*1e6);
fprintf('  PML outer size: h_out = %.2f μm\n', h_out*1e6);
fprintf('  PML grading: α = %.1f, β = %.1f\n', alpha, beta);

% Generate geometry
fprintf('\nGenerating geometry...\n');
try
    [node, elem, regionID, bdFlag, meta] = makeSphericalAirAndPML(...
        'L', L, 'r_SBC', r_SBC, 'r_PML', r_PML, ...
        'h_air', h_air, 'h_SiC', h_SiC, 'h_out', h_out, ...
        'alpha', alpha, 'beta', beta, 'visualize', true);
    
    fprintf('✓ Geometry generation successful\n');
catch ME
    fprintf('✗ Geometry generation failed: %s\n', ME.message);
    return;
end

% Display results
fprintf('\nGeometry Results:\n');
fprintf('  Total nodes: %d\n', meta.nnodes);
fprintf('  Total elements: %d\n', meta.nelems);
fprintf('  SiC elements: %d (%.1f%%)\n', meta.nSiC, 100*meta.nSiC/meta.nelems);
fprintf('  Air elements: %d (%.1f%%)\n', meta.nAir, 100*meta.nAir/meta.nelems);
fprintf('  PML elements: %d (%.1f%%)\n', meta.nPML, 100*meta.nPML/meta.nelems);
fprintf('  SBC boundary faces: %d\n', meta.nSBC);

% Validate geometry
fprintf('\nGeometry Validation:\n');

% Check region nesting
if meta.nSiC > 0 && meta.nAir > 0 && meta.nPML > 0
    fprintf('✓ All three regions present\n');
else
    fprintf('✗ Missing regions\n');
end

% Check SBC boundary
if meta.nSBC > 0
    fprintf('✓ SBC boundary faces identified\n');
else
    fprintf('⚠ No SBC boundary faces found\n');
end

% Check mesh quality
avg_element_volume = (4/3*pi*r_PML^3) / meta.nelems;
avg_element_size = avg_element_volume^(1/3);
fprintf('  Average element size: %.2f μm\n', avg_element_size*1e6);
fprintf('  Target air mesh size: %.2f μm\n', h_air*1e6);

if avg_element_size <= 2*h_air
    fprintf('✓ Mesh density adequate\n');
else
    fprintf('⚠ Mesh may be too coarse\n');
end

% Save results
fprintf('\nSaving results...\n');
save('phase1_geometry_results.mat', 'node', 'elem', 'regionID', 'bdFlag', 'meta', ...
     'L', 'r_SBC', 'r_PML', 'lambda', 'h_air', 'h_SiC', 'h_out', 'alpha', 'beta');

fprintf('✓ Results saved to phase1_geometry_results.mat\n');

fprintf('\n=== Phase 1 Complete ===\n');
fprintf('Ready for Phase 2: Scattering Boundary Condition\n');
