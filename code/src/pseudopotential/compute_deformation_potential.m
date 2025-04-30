function [D] = compute_deformation_potential(n, k, eps, deps, par)

D = zeros(3,3);

for i = 1 : 3
  ei    = zeros(3,1);
  ei(i) = 1;
  for j = 1 : 3
    ej    = zeros(3,1);
    ej(j) = 1;

    % central finite difference approx
    % add strain
      EPS = eps + 0.25*(ei*ej' + ej*ei') * deps;
      H = Hamiltonian(k, EPS, par);
      [~,E] = eig(H);
      E  = diag(E) * par.energy_scale;
      Ep = E(n);

    % subtract strain
      EPS = eps - 0.25*(ei*ej' + ej*ei') * deps;
      H = Hamiltonian(k, EPS, par);
      [~,E] = eig(H);
      E  = diag(E) * par.energy_scale;
      Em = E(n);

    % deformation
      D(i,j) = (Ep-Em)/deps;


  end
end