function [cent, Enorm] = barycenter_Enorm(node, elem, u)
% |E_h| at tetra barycenters for ND1 edge elements (no reliance on elem2edgeSign)

NT   = size(elem,1);
cent = (node(elem(:,1),:) + node(elem(:,2),:) + ...
        node(elem(:,3),:) + node(elem(:,4),:)) / 4;

[elem2edge, ~] = dof3edge(elem);          % element -> 6 global edge IDs
[DLambda, ~]   = gradbasis3(node, elem);  % NT x 3 x 4 (∇L1..∇L4)

pairs = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];   % local edges (i<j)
Enorm = zeros(NT,1);

for K = 1:NT
    grads = squeeze(DLambda(K,:,:));      % 3x4
    % ND1 basis at barycenter: N_ij = (1/4)(∇L_j - ∇L_i)
    Nloc  = zeros(3,6);
    for e = 1:6
        i = pairs(e,1);  j = pairs(e,2);
        Nloc(:,e) = 0.25*(grads(:,j) - grads(:,i));
    end
    eg  = elem2edge(K,:);                 % 1x6 global edge indices

    % Local->global edge orientation sign via node IDs (min->max is +1)
    sgn = zeros(6,1);
    for e = 1:6
        vi = elem(K, pairs(e,1)); 
        vj = elem(K, pairs(e,2));
        sgn(e) = 1; if vi > vj, sgn(e) = -1; end
    end

    ueg  = u(eg);               % 1x6 or 6x1 edge DOFs
    coef = ueg(:) .* sgn;       % 6x1
    Eh   = Nloc * coef;         % 3x1 E at barycenter (complex OK)
    Enorm(K) = norm(Eh);        % |E| (phasor magnitude)
end
end
