%QUICK_TEST Quick test to verify the setup works
%
% This script tests if all the necessary functions are available
% and provides a simple geometry generation test.
%
% Author: iFEM-Maxwell-PML Project
% Date: 2025

clear; clc; close all;

fprintf('=== Quick Setup Test ===\n\n');

% Setup paths
fprintf('Setting up paths...\n');
project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(project_root);
addpath(fullfile(project_root, 'src'));
addpath(fullfile(project_root, 'src', 'geometry'));
addpath(fullfile(project_root, 'src', 'physics'));
addpath(fullfile(project_root, 'src', 'post'));
addpath(fullfile(project_root, 'src', 'examples'));

% Check if main function exists
fprintf('Checking function availability...\n');
if exist('makeSphericalAirAndPML', 'file') == 2
    fprintf('✓ makeSphericalAirAndPML function found\n');
else
    fprintf('✗ makeSphericalAirAndPML function not found\n');
    fprintf('  Current path:\n');
    path_info = path;
    fprintf('  %s\n', path_info);
    return;
end

% Check dependencies
fprintf('Checking dependencies...\n');
if exist('distmesh3d', 'file') == 2
    fprintf('✓ distmesh3d available\n');
else
    fprintf('⚠ distmesh3d not found - will use fallback\n');
end

if exist('uniformrefine3', 'file') == 2
    fprintf('✓ uniformrefine3 available\n');
else
    fprintf('⚠ uniformrefine3 not found - refinement may not work\n');
end

if exist('findboundary3', 'file') == 2
    fprintf('✓ findboundary3 available\n');
else
    fprintf('⚠ findboundary3 not found - boundary detection may not work\n');
end

% Try a simple function call
fprintf('\nTesting function call...\n');
try
    % Simple test with minimal parameters
    [node, elem, regionID, bdFlag, meta] = makeSphericalAirAndPML(...
        'L', 1e-6, 'r_SBC', 6e-6, 'r_PML', 8e-6, ...
        'h_air', 0.5e-6, 'h_SiC', 0.25e-6, 'h_out', 1.0e-6, ...
        'visualize', false);
    
    fprintf('✓ Function call successful!\n');
    fprintf('  Generated %d nodes and %d elements\n', size(node,1), size(elem,1));
    fprintf('  SiC: %d, Air: %d, PML: %d elements\n', ...
            sum(regionID==1), sum(regionID==2), sum(regionID==3));
    fprintf('  SBC boundary faces: %d\n', sum(bdFlag==201));
    
catch ME
    fprintf('✗ Function call failed: %s\n', ME.message);
    fprintf('  Error location: %s (line %d)\n', ME.stack(1).file, ME.stack(1).line);
    return;
end

fprintf('\n=== Setup Test Complete ===\n');
fprintf('You can now run the full test scripts!\n');
