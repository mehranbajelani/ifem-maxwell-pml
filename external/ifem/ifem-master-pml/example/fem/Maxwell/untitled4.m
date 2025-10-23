% ==== geometry / mesh (use your existing cube mesh maker) ====
[node,elem] = cubemesh([-1,1,-1,1,-1,1],1);
[edge, elem2edge] = dof3edge(elem);    % global edges

% ==== constants (SI) ====
eps0  = 8.854187817e-12;   % F/m
mu0   = 4*pi*1e-7;         % H/m
er    = 1.0;               % relative permittivity
mur   = 1.0;
f     = 5e6;               % 5 MHz (coercive regime)
omega = 2*pi*f;
eta   = 1e-3;              % tiny loss (0.1%) to keep things strictly coercive

% ==== PDE struct ====
pde.mu      = mu0*mur;                 % physical mu
pde.epsilon = eps0*er*(1 - 1i*eta);    % physical epsilon (with small loss)
pde.omega   = omega;

J0 = 1;  % A/m^2
pde.J = @(P) [ J0*sin(pi*P(:,1)).*sin(pi*P(:,2)).*sin(pi*P(:,3)), ...
               0*P(:,1), 0*P(:,1) ];

% PEC (Dirichlet) on all faces
bdFlag = setboundary3(node, elem, 'Dirichlet');
pde.g_D = [];   % PEC => zero tangential E (enforced strongly by masking)

% ==== solver options ====
option.solver    = 'amg';   % AMG/HX branch
option.outsolver = 'cg';
option.tol       = 1e-8;
option.maxit     = 2000;

% ==== solve ====
[u, edge, eqn, info] = Maxwell1(node, elem, bdFlag, pde, option);
fprintf('iters=%d  stop=%.2e\n', info.itStep, info.stopErr);

% ==== global norms to compare with COMSOL ====
zeroE = @(P) [0*P(:,1), 0*P(:,1), 0*P(:,1)];
E_L2  = getL2error3ND1(node, elem, zeroE, u);  % sqrt(∫ |E|^2)
[~, volT] = gradbasis3(node, elem);
V = sum(volT);
E_rms = E_L2 / sqrt(V);

fprintf('||E||_L2  = %.6e  V/m*sqrt(m^3)\n', E_L2);
fprintf('RMS(|E|)  = %.6e  V/m\n', E_rms);
