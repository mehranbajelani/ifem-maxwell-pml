clear; clc;

% --- Physics & units ---
f   = 10e9;                             % Hz (change freely)
w   = 2*pi*f;
[mu0,eps0,c0,eta0] = physconst_EM();
k0  = w*sqrt(mu0*eps0);

% --- Problem params ---
L      = 1e-2;                          % cube side [m]
epsr   = 4.0;                           % uniform dielectric
mur    = 1.0;
sigma  = 0.0;                           % S/m (set small >0 if near resonance)
J0     = 1.0;                           % A/m^2 amplitude for impressed J

% --- Mesh ---
nsub   = 4;                              % subdivision level (refine later)
[node, elem, bdFaceId] = make_cube_mesh(L, nsub);   % tags 6 faces

% --- Edge DoFs & PEC edges ---
[edge, elem2edge] = dof3edge(elem);
isPECedge = find_pec_edges(node, edge, bdFaceId);   % logical Nedges x 1

% --- Materials per element ---
NT = size(elem,1);
epsr_elem = epsr*ones(NT,1);
mur_elem  = mur*ones(NT,1);
sigma_elem= sigma*ones(NT,1);

% --- Source (impressed J) ---
Jfun = @(x) gaussian_Jx(x, L/20, J0);   % centered compact x-directed current

% --- Assembly ---
K       = assemble_curlcurl(node, elem, elem2edge, mur_elem);
M_eps   = assemble_mass_eps(node, elem, elem2edge, epsr_elem);
M_I     = assemble_mass_I(node, elem, elem2edge);   % optional identity mass
b       = 1i*w*mu0 * assemble_rhs_J(node, elem, elem2edge, Jfun);

A = K - (k0^2)*M_eps - 1i*w*mu0*(sigma_elem(1))*M_I; % uniform sigma for now

% --- Apply PEC (Dirichlet on boundary edges) ---
[iiFree, Ared, bred] = apply_dirichlet_edges(A, b, isPECedge);

% --- Solve ---
Efree = Ared \ bred;

% --- Recover full vector (PEC edges = 0) ---
NE    = size(edge,1);
Eall  = zeros(NE,1); Eall(iiFree) = Efree;

% --- Quick checks ---
fprintf('||E||2 = %.3e\n', norm(Eall));
fprintf('Max |E| on PEC edges = %.3e\n', max(abs(Eall(isPECedge))));
