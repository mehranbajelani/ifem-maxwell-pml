%% CUBEMAXWELL11 solves Maxwell type equations in a cube using linear element.
%

close all;

%% Defacult setting
[node,elem] = cubemesh([-1,1,-1,1,-1,1],1);
pde = Maxwelldata2;
% pde = planewavedata1;
% bdFlag = setboundary3(node,elem,'Neumann');
bdFlag = setboundary3(node,elem,'Dirichlet');
% option.solver = 'mg';
option.solver = 'amg';
% showmesh3(node,elem); view([130,28]);

%% Parameters
maxIt = 4; 
N = zeros(maxIt,1); 
h = zeros(maxIt,1);
energyErr = zeros(maxIt,1);
uIuhErr = zeros(maxIt,1);
L2Err = zeros(maxIt,1);

%% Finite Element Method        
for k = 1:maxIt
    % refine grid    
    [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);
    % solve the equation
[u,edge,eqn,info] = Maxwell1(node,elem,bdFlag,pde,option);


[cent, Enorm] = barycenter_Enorm(node, elem, u);
scatter3(cent(:,1),cent(:,2),cent(:,3),12,Enorm,'filled'); colorbar
title('|E| at tetra barycenters')
title ('[E] at tetra') 
zero = @(p) [0*p(:,1), 0*p(:,1), 0*p(:,1)];
E_L2 = getL2error3ND1(node, elem, zero, u);   % equals ||E_h||_{L2(Ω)}
disp(E_L2)



    % compute error
    energyErr(k) = getHcurlerror3ND1(node,elem,pde.curlu,real(u));
    L2Err(k) = getL2error3ND1(node,elem,pde.exactu,real(u));
    uI = edgeinterpolate1(pde.exactu,node,edge);
    uIuhErr(k) = sqrt((u-uI)'*eqn.A*(u-uI));    
    N(k) = length(u);
    h(k) = 1./(size(node,1)^(1/3)-1);        
end
%% Plot convergence rates
figure;
showrateh3(h,energyErr,1,'k-+','|| curl (u-u_h) ||',...
           h,uIuhErr,1,'r-+','|| curl (u_I-u_h) ||',...
           h,L2Err,1,'b-+','|| u-u_h||');

E_L2_exact = sqrt(8 + 4*sin(2));
% Compare to your computed E_L2:
relerr = abs(E_L2 - E_L2_exact)/E_L2_exact;

% After solving and calling:
[cent, Enorm] = barycenter_Enorm(node, elem, u);

% Per-tet volumes (for averages/integrals)
v1  = node(elem(:,2),:) - node(elem(:,1),:);
v2  = node(elem(:,3),:) - node(elem(:,1),:);
v3  = node(elem(:,4),:) - node(elem(:,1),:);
vol = abs(dot(v1, cross(v2, v3, 2), 2)) / 6;   % NT×1

% Summary stats of |E|
E_min  = min(Enorm);
E_max  = max(Enorm);
E_mean = sum(Enorm .* vol) / sum(vol);                 % volume-weighted mean |E|
E_rms  = sqrt( sum((Enorm.^2).*vol) / sum(vol) );      % RMS |E|
E_L2b  = sqrt( sum((Enorm.^2).*vol) );                 % ≈ ||E||_{L2(Ω)} (barycenter quad)

disp([E_min, E_mean, E_rms, E_max, E_L2b])
