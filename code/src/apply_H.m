function [H_psi] = apply_H(psi, x, par)

  % action of kinetic energy operator
    prefactor = + par.const.hbar*par.const.hbar/(2*par.mass_l);
    T_psi = prefactor * ifft(par.k.^2 .* fft(psi));

  %%%%%%%%%%%%%%  
  % potential energy operator
    U_x = potential_modification(x, par);
    U   = par.U_QW + par.U_F + U_x;

  % action of U
    U_psi = U.*psi;

  % action of Hamiltonian    
    H_psi = (T_psi + U_psi)/par.energy_scale;


end