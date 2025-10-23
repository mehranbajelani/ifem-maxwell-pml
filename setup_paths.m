%SETUP_PATHS Add necessary paths for iFEM-Maxwell-PML project
%
% This script adds the required directories to the MATLAB path
% so that all functions can be found and executed properly.
%
% Author: iFEM-Maxwell-PML Project
% Date: 2025

fprintf('Setting up paths for iFEM-Maxwell-PML project...\n');

% Get the current directory (should be the project root)
project_root = pwd;
fprintf('Project root: %s\n', project_root);

% Add necessary directories to path
addpath(project_root);
addpath(fullfile(project_root, 'src'));
addpath(fullfile(project_root, 'src', 'geometry'));
addpath(fullfile(project_root, 'src', 'physics'));
addpath(fullfile(project_root, 'src', 'post'));
addpath(fullfile(project_root, 'src', 'examples'));

% Check if iFEM is available (for uniformrefine3, findboundary3, etc.)
if exist('uniformrefine3', 'file') == 2
    fprintf('✓ iFEM functions available\n');
else
    fprintf('⚠ iFEM functions not found - some features may not work\n');
    fprintf('  Consider adding iFEM to your path or installing it\n');
end

% Check if distmesh3d is available
if exist('distmesh3d', 'file') == 2
    fprintf('✓ distmesh3d available\n');
else
    fprintf('⚠ distmesh3d not found - will use fallback mesh generation\n');
end

% Check if our main function is available
if exist('makeSphericalAirAndPML', 'file') == 2
    fprintf('✓ makeSphericalAirAndPML function found\n');
else
    fprintf('✗ makeSphericalAirAndPML function not found\n');
    fprintf('  Check that src/geometry/ is in the path\n');
end

fprintf('Path setup complete!\n');
fprintf('You can now run: cube_sphAir_sphPML_driver.m or test_geometry.m\n');
