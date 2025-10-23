%% CUBEMAXWELLSADDLE solves Maxwell type equations in a cube using linear order element.
% This is a special case of div u = g being nozero.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

clear; close all;

%% Defacult setting
[node,elem] = cubemesh([-1,1,-1,1,-1,1],1);
%%
pde.J = @(p) [sin(p(:,1)).*cos(p(:,2)).*sin(p(:,3)), ...
              cos(p(:,1)).*sin(p(:,2)).*sin(p(:,3)), ...
              2*cos(p(:,1)).*cos(p(:,2)).*cos(p(:,3))];
pde.exactu = @(p)[0*p(:,1), 0*p(:,2), ...
    cos(p(:,1)).*cos(p(:,2)).*cos(p(:,3))];
pde.g_D = pde.exactu;
pde.curlu = @(p) [-cos(p(:,1)).*sin(p(:,2)).*cos(p(:,3)), ...
              sin(p(:,1)).*cos(p(:,2)).*cos(p(:,3)), 0*p(:,3)];
pde.g = @(p) -cos(p(:,1)).*cos(p(:,2)).*sin(p(:,3));
pde.mu = 1;

%%
bdFlag = setboundary3(node,elem,'Dirichlet');
option.printlevel = 0;
% option.solver = 'direct';
option.solver = 'mg';
option.solver = 'diag';

%% Parameters
maxIt = 3; 
N = zeros(maxIt,1); 
h = zeros(maxIt,1);
energyErr = zeros(maxIt,1);
L2Err = zeros(maxIt,1);
uIuhErr = zeros(maxIt,1);

%% Finite Element Method        
for k = 1:maxIt   
    [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);     
    [soln,eqn,info] = Maxwell1saddle(node,elem,bdFlag,pde,option);  
    u = soln.u;
    fprintf('\n\n # of DoFs = %d \n',length(u));
    % compute error
    uI = edgeinterpolate1(pde.exactu,node,eqn.edge);
    energyErr(k) = getHcurlerror3ND1(node,elem,pde.curlu,u);
    L2Err(k) = getL2error3ND1(node,elem,pde.exactu,u);
    uIuhErr(k) = sqrt((u-uI)'*(eqn.A)*(u-uI));        
%     L2Err(k) = sqrt((u-uI)'*(eqn.M)*(u-uI));        
    fprintf('\n ||curl(u-u_h)|| is %g \n',energyErr(k))
    N(k) = length(u);
    h(k) = 1./(size(node,1)^(1/3)-1);   
end

%% Plot convergence rates
figure(1);
showrateh3(h,energyErr,1,'k-+','|| curl (u-u_h) ||',...
           h,uIuhErr,1,'r-+','|| curl (u_I-u_h) ||',...
           h,L2Err,1,'b-+','|| u-u_h||');