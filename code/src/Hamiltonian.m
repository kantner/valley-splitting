function H = Hamiltonian(k_vec, par)
% Hamiltonian describing Si bandstructure using non-local empirical pseudopotentials
% References:
%   Fischetti and Laux, J. Appl. Phys. 80, 2234-2252 (1996)
%   Fischetti and Higman, Theory and calculation of the deformation potential electron-phonon scattering rates in semiconductors. In: Monte Carlo Device Simulation, Hess, Karl (Ed.), p. 123-160, Springer (1991)


% imports
  N_G     = par.N_G;
  G_list  = par.G_list;
  tau     = par.tau;
  Omega_a = par.Omega_a;

% constants
  hbar = par.const.hbar;  
  m0   = par.const.m0;

% Fermi vector and Fermi energy
  KF = par.KF;
  EF = (hbar*KF)^2/(2*m0);

% nonlocal 
  non_local_switch = 1; % 0 = off | 1 = on
  V_nloc = 0;

% scaling of energy
  energy_scale = par.energy_scale;
  
%%%%%%%%%%%%%%%%%%%%%  
% construct Hamiltonian
  H = zeros(N_G,N_G);
  
%%%%%%%%%%%%%%%%%%%%%
% kinetic energy
  prefactor = hbar*hbar/(2*m0);
  for i = 1 : N_G
    H(i,i) = prefactor * (k_vec + G_list(:,i))' * (k_vec + G_list(:,i)) * 1/energy_scale;
  end


%%%%%%%%%%%%%%%%%%%%%
% potential energy
  for i = 1 : N_G

    %%%%%%%%%%%
    % diagonal terms
    % structure factor S = cos(0.5*tau*(G-G')) is S(0) = 1 for G=G'
      S = 1;
    
    % local potential at G-G'=0 (has to be evaluated only once)
      if i == 1 
        V0 = V_atomic(0, par);
      end
      V_loc = S * V0  * 1/energy_scale;

    % non-local potential
      if non_local_switch == 1
        Ki      = k_vec + G_list(:,i);
        norm_Ki = sqrt(Ki'*Ki);
        E1      = prefactor * norm_Ki * norm_Ki;
        A       = par.alpha0 + par.beta0 * (E1 - EF);
        F       = function_F(norm_Ki,norm_Ki, par.R0);

        V_nloc =  4*pi/Omega_a * A * S * F  * 1/energy_scale;
      end

    % add to kinetic energy  
      H(i,i) = H(i,i) + V_loc + V_nloc;

    %%%%%%%%%%%%%%%%%%%%%% 
    % off-diagonal
    for j = i+1 : N_G

    % difference vector  
      Gij = G_list(:,i) - G_list(:,j);

    % structure factor
      S = cos(0.5*tau'*Gij);
      
    % local potential  
      norm_Gij = sqrt(Gij'*Gij);
      V_a      = V_atomic(norm_Gij, par);
      V_loc    = S * V_a  * 1/energy_scale;

    % non-local potential
      if non_local_switch == 1
        Ki = k_vec + G_list(:,i);
        Kj = k_vec + G_list(:,j);

        norm_Ki2 = Ki'*Ki;
        norm_Kj2 = Kj'*Kj;

        norm_Ki = sqrt(norm_Ki2);
        norm_Kj = sqrt(norm_Kj2);

        E1 = prefactor * norm_Ki2;
        E2 = prefactor * norm_Kj2;
        A  = par.alpha0 + par.beta0 * (sqrt(E1*E2)/par.units.Ry - EF/par.units.Ry) * par.units.Ry;
        F = function_F(norm_Ki,norm_Kj, par.R0);

        V_nloc = 4*pi/Omega_a * A * S * F * 1/energy_scale;
      end

    % add to Hamiltonian
      H(i,j) = V_loc + V_nloc;

      if isnan(H(i,j))
        pause
      end

    % add conjugate term  
      H(j,i)   = H(i,j)';
      

    end

  end

end