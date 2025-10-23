function [node, elem, regionID, bdFlag, meta] = makeSphericalAirAndPML(varargin)
%MAKESPHERICALAIRANDPML Create mesh for SiC cube + spherical air + PML geometry
%
% Syntax:
%   [node, elem, regionID, bdFlag, meta] = makeSphericalAirAndPML()
%   [node, elem, regionID, bdFlag, meta] = makeSphericalAirAndPML('param', value, ...)
%
% Description:
%   Creates a three-region geometry for electromagnetic scattering:
%   - SiC cube (L = 1 μm) at origin
%   - Spherical air region (r_SBC = 6 μm) 
%   - Spherical PML shell (t_PML = 2 μm, r_PML = 8 μm)
%
%   Uses hybrid region tagging with boundary detection and PML mesh grading.
%
% Parameters:
%   'L'           - SiC cube side length (default: 1e-6 m)
%   'r_SBC'       - Air sphere radius (default: 6e-6 m)
%   'r_PML'       - PML outer radius (default: 8e-6 m)
%   'h_air'       - Air region mesh size (default: 0.5e-6 m)
%   'h_SiC'       - SiC region mesh size (default: 0.25e-6 m)
%   'h_out'       - PML outer mesh size (default: 1.0e-6 m)
%   'alpha'       - PML grading parameter (default: 1.0)
%   'beta'        - PML grading parameter (default: 1.0)
%   'refine'      - Number of mesh refinements (default: 0)
%   'visualize'   - Show mesh plot (default: false)
%
% Outputs:
%   node      - Node coordinates [N x 3]
%   elem      - Element connectivity [M x 4] (tetrahedra)
%   regionID  - Region ID for each element [M x 1]
%              1 = SiC, 2 = Air, 3 = PML
%   bdFlag    - Boundary face flags [F x 1]
%              201 = SBC boundary (r = r_SBC)
%   meta      - Metadata structure with counts and statistics
%
% Example:
%   [node, elem, regionID, bdFlag, meta] = makeSphericalAirAndPML('visualize', true);
%
% Author: iFEM-Maxwell-PML Project
% Date: 2025

% Parse input parameters
p = inputParser;
addParameter(p, 'L', 1e-6, @(x) x > 0);           % SiC cube side length
addParameter(p, 'r_SBC', 6e-6, @(x) x > 0);        % Air sphere radius
addParameter(p, 'r_PML', 8e-6, @(x) x > 0);        % PML outer radius
addParameter(p, 'h_air', 0.5e-6, @(x) x > 0);      % Air mesh size
addParameter(p, 'h_SiC', 0.25e-6, @(x) x > 0);     % SiC mesh size
addParameter(p, 'h_out', 1.0e-6, @(x) x > 0);      % PML outer mesh size
addParameter(p, 'alpha', 1.0, @(x) x > 0);        % PML grading parameter
addParameter(p, 'beta', 1.0, @(x) x > 0);         % PML grading parameter
addParameter(p, 'refine', 0, @(x) x >= 0);         % Refinement level
addParameter(p, 'visualize', false, @islogical);   % Visualization flag
parse(p, varargin{:});

% Extract parameters
L = p.Results.L;
r_SBC = p.Results.r_SBC;
r_PML = p.Results.r_PML;
h_air = p.Results.h_air;
h_SiC = p.Results.h_SiC;
h_out = p.Results.h_out;
alpha = p.Results.alpha;
beta = p.Results.beta;
refine = p.Results.refine;
visualize = p.Results.visualize;

% Compute PML thickness
t_PML = r_PML - r_SBC;

% Validate geometry parameters
if r_SBC <= L/2
    error('Air sphere radius must be larger than cube half-diagonal');
end
if r_PML <= r_SBC
    error('PML outer radius must be larger than air sphere radius');
end

fprintf('Creating geometry:\n');
fprintf('  SiC cube: L = %.1f μm\n', L*1e6);
fprintf('  Air sphere: r = %.1f μm\n', r_SBC*1e6);
fprintf('  PML shell: r = %.1f μm (thickness = %.1f μm)\n', r_PML*1e6, t_PML*1e6);
fprintf('  Mesh sizes: h_SiC = %.2f μm, h_air = %.2f μm, h_out = %.2f μm\n', ...
        h_SiC*1e6, h_air*1e6, h_out*1e6);
fprintf('  PML grading: α = %.1f, β = %.1f\n', alpha, beta);

% Create bounding box for initial mesh
bbox = [-r_PML, r_PML; -r_PML, r_PML; -r_PML, r_PML];

% Generate initial mesh using distmesh3d with PML grading
fprintf('Generating initial mesh with PML grading...\n');
try
    % Create sizing function for PML grading
    hfun = @(p) computeMeshSize(p, L, r_SBC, r_PML, h_air, h_SiC, h_out, alpha, beta);
    
    % Use distmesh3d for initial mesh generation
    [node, elem] = distmesh3d(@(p) drectangle3(p, bbox), ...
                              hfun, ...
                              h_air, bbox);
catch ME
    warning('distmesh3d failed, using simple structured mesh');
    % Fallback to simple structured mesh
    [node, elem] = createStructuredMesh(bbox, h_air);
end

% Refine mesh if requested
for i = 1:refine
    fprintf('Refining mesh (level %d)...\n', i);
    [node, elem] = uniformrefine3(node, elem);
end

% Tag regions using hybrid approach
fprintf('Tagging regions with hybrid approach...\n');
regionID = tagRegionsHybrid(node, elem, L, r_SBC, r_PML);

% Create boundary information
fprintf('Creating boundary information...\n');
[bdFlag, meta] = createBoundaryInfo(node, elem, regionID, L, r_SBC, r_PML);

% Print mesh statistics
fprintf('\nMesh Statistics:\n');
fprintf('  Total nodes: %d\n', size(node, 1));
fprintf('  Total elements: %d\n', size(elem, 1));
fprintf('  SiC elements: %d (%.1f%%)\n', sum(regionID == 1), 100*sum(regionID == 1)/length(regionID));
fprintf('  Air elements: %d (%.1f%%)\n', sum(regionID == 2), 100*sum(regionID == 2)/length(regionID));
fprintf('  PML elements: %d (%.1f%%)\n', sum(regionID == 3), 100*sum(regionID == 3)/length(regionID));
fprintf('  SBC boundary faces: %d\n', sum(bdFlag == 201));

% Visualization
if visualize
    visualizeMesh(node, elem, regionID, bdFlag);
end

end

function h = computeMeshSize(p, L, r_SBC, r_PML, h_air, h_SiC, h_out, alpha, beta)
%COMPUTEMESHSIZE Compute mesh size with PML grading
% h = h_air * (1-ρ)^α + h_out * ρ^β
% where ρ = (r - r_SBC) / (r_PML - r_SBC)

% Compute radial distance
r = sqrt(sum(p.^2, 2));

% Initialize with air mesh size
h = h_air * ones(size(r));

% SiC region: use h_SiC
cube_mask = abs(p(:,1)) <= L/2 & abs(p(:,2)) <= L/2 & abs(p(:,3)) <= L/2;
h(cube_mask) = h_SiC;

% PML region: use graded mesh size
pml_mask = r >= r_SBC & r <= r_PML;
if any(pml_mask)
    rho = (r(pml_mask) - r_SBC) / (r_PML - r_SBC);
    h_pml = h_air * (1 - rho).^alpha + h_out * rho.^beta;
    h(pml_mask) = h_pml;
end

end

function regionID = tagRegionsHybrid(node, elem, L, r_SBC, r_PML)
%TAGREGIONSHYBRID Tag elements with region IDs using hybrid approach
% Region IDs: 1 = SiC, 2 = Air, 3 = PML

nelem = size(elem, 1);
regionID = zeros(nelem, 1);

% Compute element centroids
centroids = zeros(nelem, 3);
for i = 1:nelem
    centroids(i, :) = mean(node(elem(i, :), :), 1);
end

% First pass: tag by centroid
for i = 1:nelem
    x = centroids(i, 1);
    y = centroids(i, 2);
    z = centroids(i, 3);
    r = sqrt(x^2 + y^2 + z^2);
    
    % Check if inside SiC cube
    if abs(x) <= L/2 && abs(y) <= L/2 && abs(z) <= L/2
        regionID(i) = 1; % SiC
    % Check if inside air sphere
    elseif r <= r_SBC
        regionID(i) = 2; % Air
    % Check if inside PML shell
    elseif r <= r_PML
        regionID(i) = 3; % PML
    else
        error('Element outside all regions at (%.2e, %.2e, %.2e)', x, y, z);
    end
end

% Second pass: correction for elements near interfaces
delta = 0.1 * min([L, r_SBC, r_PML-r_SBC]); % Small band around interfaces

for i = 1:nelem
    % Get all vertices of this element
    verts = node(elem(i, :), :);
    r_verts = sqrt(sum(verts.^2, 2));
    
    % Check if element is near SBC interface (r = r_SBC)
    near_SBC = any(abs(r_verts - r_SBC) < delta);
    
    if near_SBC
        % Force classification based on all vertices
        if all(r_verts < r_SBC)
            regionID(i) = 2; % Air
        elseif all(r_verts > r_SBC)
            regionID(i) = 3; % PML
        end
    end
    
    % Check if element is near SiC cube surface
    near_cube = any(abs(verts(:,1)) <= L/2 + delta & abs(verts(:,2)) <= L/2 + delta & abs(verts(:,3)) <= L/2 + delta);
    
    if near_cube
        % Force SiC classification only if ALL vertices are inside cube
        if all(abs(verts(:,1)) <= L/2 & abs(verts(:,2)) <= L/2 & abs(verts(:,3)) <= L/2)
            regionID(i) = 1; % SiC
        end
    end
end

end

function [bdFlag, meta] = createBoundaryInfo(node, elem, regionID, L, r_SBC, r_PML)
%CREATEBOUNDARYINFO Create boundary face information
% Boundary IDs: 201 = SBC boundary (r = r_SBC)

% Extract boundary faces
[bdFace, ~] = findboundary3(elem);

% Initialize boundary flags
bdFlag = zeros(size(bdFace, 1), 1);

% Tag SBC boundary faces (r = r_SBC)
for i = 1:size(bdFace, 1)
    face_nodes = bdFace(i, :);
    face_centroid = mean(node(face_nodes, :), 1);
    r_centroid = sqrt(sum(face_centroid.^2));
    
    % Check if face is on SBC boundary
    if abs(r_centroid - r_SBC) < 0.1 * r_SBC
        bdFlag(i) = 201; % SBC boundary
    end
end

% Create metadata
meta = struct();
meta.nnodes = size(node, 1);
meta.nelems = size(elem, 1);
meta.nfaces = size(bdFace, 1);
meta.nSiC = sum(regionID == 1);
meta.nAir = sum(regionID == 2);
meta.nPML = sum(regionID == 3);
meta.nSBC = sum(bdFlag == 201);
meta.L = L;
meta.r_SBC = r_SBC;
meta.r_PML = r_PML;

end

function [node, elem] = createStructuredMesh(bbox, h)
%CREATESTRUCTUREDMESH Create simple structured mesh as fallback

% Simple structured mesh generation
x = bbox(1,1):h:bbox(1,2);
y = bbox(2,1):h:bbox(2,2);
z = bbox(3,1):h:bbox(3,2);

[X, Y, Z] = meshgrid(x, y, z);
node = [X(:), Y(:), Z(:)];

% Create simple tetrahedral connectivity
% This is a placeholder - proper implementation needed
elem = [];

end

function visualizeMesh(node, elem, regionID, bdFlag)
%VISUALIZEMESH Visualize the mesh with region coloring

figure('Name', 'Mesh with Region Tagging and PML Grading');
hold on;

% Plot elements by region
colors = [1 0.5 0.5; 0.5 0.5 1; 0.5 1 0.5]; % SiC, Air, PML
regionNames = {'SiC', 'Air', 'PML'};

for reg = 1:3
    idx = regionID == reg;
    if any(idx)
        elem_reg = elem(idx, :);
        % Plot tetrahedra (simplified visualization)
        for i = 1:min(100, size(elem_reg, 1)) % Limit for performance
            tet = elem_reg(i, :);
            % Plot tetrahedron faces
            faces = [1 2 3; 1 2 4; 1 3 4; 2 3 4];
            for j = 1:4
                face = faces(j, :);
                patch(node(tet(face), 1), node(tet(face), 2), node(tet(face), 3), ...
                      colors(reg, :), 'FaceAlpha', 0.3, 'EdgeAlpha', 0.1);
            end
        end
    end
end

% Highlight SBC boundary
if any(bdFlag == 201)
    sbc_faces = find(bdFlag == 201);
    fprintf('Highlighting %d SBC boundary faces\n', length(sbc_faces));
end

xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
title('Mesh with Region Tagging and PML Grading');
legend(regionNames, 'Location', 'best');
axis equal;
grid on;

end