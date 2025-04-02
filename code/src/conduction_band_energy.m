function [E] = conduction_band_energy(K_vec, par)
% input is K_vec = k_vec * a0/(2*pi)
% output is energy in eV

% rescale
  k_vec = K_vec * 2*pi/par.a0;

% solve eigenvalue problem  
  H = Hamiltonian(k_vec, par);
  energy = eig(H,'vector');

% rescale to SI units        
  E = energy(par.idx_CB_low) * par.energy_scale;

% rescale to eV  
  E = E/par.units.eV;

end