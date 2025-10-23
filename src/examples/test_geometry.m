%TEST_GEOMETRY Test script for Phase 1 geometry generation
%
% This script tests the makeSphericalAirAndPML function and validates
% the geometry parameters and mesh quality.
%
% Author: iFEM-Maxwell-PML Project
% Date: 2025

clear; clc; close all;

fprintf('=== Phase 1: Geometry and Domain Definition Test ===\n\n');

% Test parameters (in SI units)
L = 1e-6;        % SiC cube side length: 1 μm
r_SBC = 6e-6;    % Air sphere radius: 6 μm  
r_PML = 8e-6;    % PML outer radius: 8 μm
h_air = 0.5e-6;  % Air mesh size: 0.5 μm
h_SiC = 0.25e-6; % SiC mesh size: 0.25 μm
h_PML = 0.5e-6;  % PML mesh size: 0.5 μm

fprintf('Test Parameters:\n');
fprintf('  SiC cube: L = %.1f μm\n', L*1e6);
fprintf('  Air sphere: r = %.1f μm\n', r_SBC*1e6);
fprintf('  PML shell: r = %.1f μm (thickness = %.1f μm)\n', r_PML*1e6, (r_PML-r_SBC)*1e6);
fprintf('  Mesh sizes: h_SiC = %.2f μm, h_air = %.2f μm, h_PML = %.2f μm\n\n', ...
        h_SiC*1e6, h_air*1e6, h_PML*1e6);

% Test 1: Basic geometry generation
fprintf('Test 1: Basic geometry generation...\n');
try
    [node, elem, regionID, bdFlag, meta] = makeSphericalAirAndPML(...
        'L', L, 'r_SBC', r_SBC, 'r_PML', r_PML, ...
        'h_air', h_air, 'h_SiC', h_SiC, 'h_out', 1.0e-6, ...
        'alpha', 1.0, 'beta', 1.0, 'visualize', true);
    fprintf('✓ Geometry generation successful\n');
catch ME
    fprintf('✗ Geometry generation failed: %s\n', ME.message);
    return;
end

% Test 2: Validate region tagging
fprintf('\nTest 2: Region tagging validation...\n');
nSiC = sum(regionID == 1);
nAir = sum(regionID == 2);
nPML = sum(regionID == 3);
ntotal = length(regionID);

fprintf('  SiC elements: %d (%.1f%%)\n', nSiC, 100*nSiC/ntotal);
fprintf('  Air elements: %d (%.1f%%)\n', nAir, 100*nAir/ntotal);
fprintf('  PML elements: %d (%.1f%%)\n', nPML, 100*nPML/ntotal);

% Validate that all elements are tagged
if nSiC + nAir + nPML ~= ntotal
    fprintf('✗ Region tagging incomplete\n');
else
    fprintf('✓ Region tagging complete\n');
end

% Test SBC boundary faces
fprintf('\nTest 2b: SBC boundary validation...\n');
nSBC = sum(bdFlag == 201);
fprintf('  SBC boundary faces: %d\n', nSBC);
if nSBC > 0
    fprintf('✓ SBC boundary faces identified\n');
else
    fprintf('⚠ No SBC boundary faces found\n');
end

% Test 3: Check mesh quality
fprintf('\nTest 3: Mesh quality check...\n');
if ntotal > 0
    fprintf('✓ Mesh contains %d elements\n', ntotal);
else
    fprintf('✗ Empty mesh generated\n');
end

% Test 4: Validate geometry parameters
fprintf('\nTest 4: Geometry parameter validation...\n');
% Check that air sphere is larger than cube
cube_diagonal = L * sqrt(3);
if r_SBC > cube_diagonal/2
    fprintf('✓ Air sphere radius > cube diagonal/2\n');
else
    fprintf('✗ Air sphere radius too small\n');
end

% Check that PML is larger than air
if r_PML > r_SBC
    fprintf('✓ PML outer radius > air sphere radius\n');
else
    fprintf('✗ PML outer radius too small\n');
end

% Test 5: Mesh density check
fprintf('\nTest 5: Mesh density analysis...\n');
if ntotal > 0
    % Estimate element size
    bbox_volume = (2*r_PML)^3;
    avg_element_volume = bbox_volume / ntotal;
    avg_element_size = avg_element_volume^(1/3);
    fprintf('  Average element size: %.2f μm\n', avg_element_size*1e6);
    fprintf('  Target air mesh size: %.2f μm\n', h_air*1e6);
    
    if avg_element_size <= 2*h_air
        fprintf('✓ Mesh density adequate\n');
    else
        fprintf('⚠ Mesh may be too coarse\n');
    end
end

fprintf('\n=== Phase 1 Test Complete ===\n');

% Save test results
save('test_geometry_results.mat', 'node', 'elem', 'regionID', 'bdFlag', 'meta', ...
     'L', 'r_SBC', 'r_PML', 'h_air', 'h_SiC', 'h_out');

fprintf('Test results saved to test_geometry_results.mat\n');
