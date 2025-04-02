function [mass] = compute_effective_mass(n, k_vec, dk, par)

  % n = band index
  Hessian = zeros(3,3);

  for i = 1 : 3
      ei    = zeros(3,1);
      ei(i) = 1;
    
    for j = i : 3
      ej    = zeros(3,1);
      ej(j) = 1;

    % perturbations  
      H  = Hamiltonian(k_vec + 0.5*dk*(ei+ej), par);
      E = eig(H,'vector');
      E1 = E(n);

      H  = Hamiltonian(k_vec + 0.5*dk*(ej-ei), par);
      E = eig(H,'vector');
      E2 = E(n);

      H  = Hamiltonian(k_vec - 0.5*dk*(ej-ei), par);
      E = eig(H,'vector');
      E3 = E(n);

      H  = Hamiltonian(k_vec - 0.5*dk*(ej+ei), par);
      E = eig(H,'vector');
      E4 = E(n);

    % compute Hessian  
      Hessian(i,j) = (E1 - E2 - E3 + E4)/(dk*dk);
      if i~=j
        Hessian(j,i) = Hessian(i,j);
      end


    end
  end

  hbar = par.const.hbar;
  m0   = par.const.m0;

% effective mass tensor (SI units)  
  mass = hbar*hbar*inv(Hessian) / par.energy_scale;

% effective mass tensor (in electron masses)  
  mass = mass/m0;




end