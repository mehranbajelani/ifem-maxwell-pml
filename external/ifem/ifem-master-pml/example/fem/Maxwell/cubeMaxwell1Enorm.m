close all;
[node,elem] = cubemesh([-1,1,-1,1,-1,1],1);

% --- constants & frequency ---
c0   = 299792458;
mu0  = 4*pi*1e-7;
eta0 = mu0*c0;
f0   = 1e9;
k0   = 2*pi*f0/c0;     % <-- wavenumber

% --- material (loss via Im(epsr), no sigma) ---
pde.mu      = 1.0;
pde.epsilon = 4.0 - 1.0i;
pde.omega   = k0;      % <-- Maxwell1 expects k0 (1/m)


% Gaussian x-directed current at center
% Physical domain current density J_e (A/m^2), x-directed Gaussian:
sigmaJ = 0.25;                   % matches COMSOL's 0.25 m in [-1,1]^3 m
J0 = 1;                          % amplitude A/m^2
Jfun = @(x) [ J0*exp(-(x(:,1).^2 + x(:,2).^2 + x(:,3).^2)/sigmaJ^2), ...
              zeros(size(x,1),1), zeros(size(x,1),1) ];

% RHS = j*omega*mu0 * J_e  = j*k0*eta0 * J_e
pde.J = @(x) 1i*k0*eta0 * Jfun(x);
% pde.f   = @(x) 1i*k0*eta0 * Jfun(x);   % <— use this in Maxwell1 if it expects 'f'
pde.g_D = @(x) zeros(size(x,1),3);     % PEC: zero tangential E

% Make Dirichlet non-empty so Maxwell1 keeps bdFlag (PEC)
pde.g_D = @(x) zeros(size(x,1),3);
pde.g_N = []; pde.g_R = [];

% PEC on all faces and pass through refinements
bdFlag = setboundary3(node,elem,'Dirichlet');

% AMG + non-CG Krylov (Maxwell1 will now respect this after your edit)
option.solver     = 'amg';
option.outsolver  = 'gmres';
option.printlevel = 1;
option.tol        = 1e-8;
option.maxiter   = 8000;                   % give it room, just for this test

maxIt = 3; N = zeros(maxIt,1); h = zeros(maxIt,1);
for k = 1:maxIt
    if k==2, option.solver='direct'; end   % quick sanity: ~8k DOFs is fine for direct

    [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);
    [u,edge,eqn] = Maxwell1(node,elem,bdFlag,pde,option);
    fprintf('#%d: dof=%d, ||b||=%.3e, ||u||=%.3e\n', k, length(u), norm(eqn.f), norm(u));
    N(k) = length(u);  h(k) = 1./(size(node,1)^(1/3)-1);
    relres = norm(eqn.bigA*u - eqn.f)/norm(eqn.f);
Ecurl  = real(u'*(eqn.A)*u);      % ~ ∥curl E∥² term
Emass  = real(u'*(eqn.M)*u);      % ~ ω²⟨εE,E⟩ term
fprintf('relres=%.2e,  ||u||=%.3e,  Ecurl=%.3e,  Emass=%.3e\n', ...
        relres, norm(u), Ecurl, Emass);

end

% --- Edge-averaged tangential |E| estimate from edge DOFs ---
% A Nedelec DOF is the line integral \int_edge E·t ds ≈ (avg E_t)*|edge|.
% So avg |E_t| on an edge ≈ |u_e| / |edge|.
NE   = size(edge,1);
elen = sqrt(sum((node(edge(:,1),:) - node(edge(:,2),:)).^2,2));   % NE x 1
u_phi = u(1:NE);             % curl-carrying block
u_psi = u(NE+1:2*NE);        % curl-free block

% edge-avg |E_t| proxy (RSS of the two edge coefficients, scaled by edge length)
E_edge = hypot(abs(u_phi)./elen, abs(u_psi)./elen);

mid = (node(edge(:,1),:) + node(edge(:,2),:))/2;                  % NE x 3


% central-box stats (where the source sits)
maskC = all(abs(mid) <= 0.25, 2);
fprintf('center box: median |E|_edge = %.3e, max = %.3e\n', ...
        median(E_edge(maskC)), max(E_edge(maskC)));

% 3D scatter of edges colored by |E| (edge-avg)
figure; scatter3(mid(:,1), mid(:,2), mid(:,3), 6, E_edge, 'filled');
axis equal tight; view(3); colorbar; title('|E_t| edge-avg');
xlabel('x'); ylabel('y'); zlabel('z');

% 2D slice at z ≈ 0 (edges whose midpoints lie on z=0 plane)
maskZ = abs(mid(:,3)) < 1e-8;
figure; scatter(mid(maskZ,1), mid(maskZ,2), 10, E_edge(maskZ), 'filled');
axis equal tight; colorbar; title('|E_t| edge-avg on z=0 slice');
xlabel('x'); ylabel('y');

if isfield(eqn,'isBdEdge')
    maxPEC = max([abs(u(eqn.isBdEdge)); abs(u(eqn.isBdEdge + NE))]);
    fprintf('max |u| on PEC edges (both blocks) = %.3e\n', maxPEC);
end
V = 8; % volume of [-1,1]^3
Erms_A = sqrt( Emass / ( (pde.omega^2) * real(pde.epsilon) * V ) ); % if M includes k0^2
Erms_B = sqrt( Emass / ( real(pde.epsilon) * V ) );                 % if M excludes k0^2
fprintf('Erms_A=%.6e  Erms_B=%.6e  (A assumes M has k0^2)\n', Erms_A, Erms_B);


%%
% --- ∫ |E| dV using the existing ND1 machinery (drop in after the solve) ---
[elem2dof,~] = dof3edge(elem);                          % DOF map (two blocks)  :contentReference[oaicite:3]{index=3}
NT   = size(elem,1);  Ndof = max(elem2dof(:));
[Dlambda,vol] = gradbasis3(node,elem);                  % grads + element volumes :contentReference[oaicite:4]{index=4}
[lambda,w]    = quadpts3(3);                            % same quad order as helper :contentReference[oaicite:5]{index=5}
locEdge = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];

acc = zeros(NT,1);
for p = 1:size(lambda,1)
    Ehp = zeros(NT,3);
    for k = 1:6
        i = locEdge(k,1);  j = locEdge(k,2);
        % φ-block and ψ-block contributions, exactly as in getL2error3ND1:  :contentReference[oaicite:6]{index=6}
        Ehp = Ehp + repmat(u(elem2dof(:,k)),1,3)         .* ( lambda(p,i)*Dlambda(:,:,j) - lambda(p,j)*Dlambda(:,:,i) );
        Ehp = Ehp + repmat(u(elem2dof(:,k)+Ndof),1,3)    .* ( lambda(p,i)*Dlambda(:,:,j) + lambda(p,j)*Dlambda(:,:,i) );
    end
    acc = acc + w(p)*sqrt(sum(abs(Ehp).^2,2));           % |E| at quadrature → sum per element
end
int_absE   = sum(acc .* vol);                            % ∫ |E| dV
avg_absE   = int_absE / sum(vol);                        % volume-average of |E|
fprintf('∫|E| dV = %.6e,   avg|E| = %.6e V/m\n', int_absE, avg_absE);

%% graphics 
% ---------- COMSOL-style |E| slice on z = z0 (robust, high-quality) ----------
z0 = 0;                % slice plane
nx = 401; ny = 401;    % increase for smoother images

% (1) build a grid that's offset half a cell so we don't sample right on faces
xlim = [min(node(:,1)) max(node(:,1))];
ylim = [min(node(:,2)) max(node(:,2))];
dx = diff(xlim)/nx;    dy = diff(ylim)/ny;
xv = linspace(xlim(1)+dx/2, xlim(2)-dx/2, nx);
yv = linspace(ylim(1)+dy/2, ylim(2)-dy/2, ny);
[xg,yg] = ndgrid(xv,yv);

% tiny z jitter avoids exact coplanar evaluations (degenerate barycentrics)
zJitter = 1e-10;
q = [xg(:), yg(:), (z0+zJitter)*ones(nx*ny,1)];

% triangulation & location
TR = triangulation(elem,node);
[ti, bc] = pointLocation(TR, q);
inside = ~isnan(ti);

% precompute element data
[elem2dof,~] = dof3edge(elem);
[Dlambda,~]  = gradbasis3(node,elem);
NE   = size(edge,1);  Ndof = 2*NE;
u_phi = u(1:NE);  u_psi = u(NE+1:2*NE);
locEdge = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];

Eabs = nan(nx*ny,1);

% (2) vectorized per-element accumulation
eList = unique(ti(inside)).';
for e = eList
    I = find(ti==e);                       % sample indices in element e
    lam = bc(I,:);                         % (#I x 4)
    G   = squeeze(Dlambda(e,:,:));         % 3x4 (columns: grad L1..L4)
    Ee  = zeros(numel(I),3);
    dofs = elem2dof(e,:);
    for k = 1:6
        a = locEdge(k,1); b = locEdge(k,2);
        Phi =  lam(:,a)*(G(:,b).') - lam(:,b)*(G(:,a).');  % (#I x 3)
        Psi =  lam(:,a)*(G(:,b).') + lam(:,b)*(G(:,a).');
        Ee  = Ee + u_phi(dofs(k))*Phi + u_psi(dofs(k))*Psi;
    end
    Eabs(I) = sqrt(sum(abs(Ee).^2,2));
end

% (3) fill any outside/NaN pixels via natural-neighbor interpolation
bad = ~inside | isnan(Eabs);
if any(bad)
    Ffill = scatteredInterpolant(xg(~bad), yg(~bad), Eabs(~bad), 'natural', 'nearest');
    Eabs(bad) = Ffill(xg(bad), yg(bad));
end

% (4) optional anti-aliasing (comment out if you don't have Image Processing Toolbox)
% Eimg = imgaussfilt(reshape(Eabs,nx,ny), 0.6);
Eimg = reshape(Eabs, nx, ny);

% (5) plot
figure;
contourf(yv, xv, Eimg, 60, 'LineColor','none');  % 60 levels, smooth shading
axis equal tight; xlabel('y'); ylabel('x'); title('|E| on z=0 slice'); colorbar;
