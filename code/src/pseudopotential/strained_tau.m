function [tau] = strained_tau(eps, par)

  Id = eye(3);

%%%%%%%%%%%%%%%%
% macroscopic strain: circumcenter
  tau_0 = (Id + eps) * par.pp.tau;

%%%%%%%%%%%%%%%%  
% internal strain: barycenter (equal distance to neighboring atoms of other fcc-sublattice)

% strained lattice vectors
  for i = 1 : 3
    a{i} = (Id + eps) * par.pp.a{i};
  end

  A = [a{1}'; a{2}'; a{3}'];
  b = zeros(3,1);
  for i = 1 : 3
    b(i) = 0.5 * a{i}'*a{i};
  end

  tau_1 = A\b;

%%%%%%%%%%%%%%%%
% convex combination
  xi  = par.pp.internal_strain;
  tau = xi * tau_1 + (1-xi) * tau_0;

end