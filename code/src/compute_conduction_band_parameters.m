function [par] = compute_conduction_band_parameters(eps, xi, par)
  % compute
  % - (reciprocal) lattice vectors
  % - conduction band minimum
  % - effective mass
  % - band gap
  % - deformation potentials
  % - Fourier domain Bloch factors
  % - bandstructure coefficents
  % as a function of strain


  % check
    if sum(sum(abs(eps - eps'))) > 0
      error('strain tensor must be symmetric')
    end

  % set input
    par.eps = eps;
    par.xi  = xi;

    eps_print = par.eps * 100;
    fprintf(1,'\n')
    fprintf(1,'Compute conduction band parameters for strain tensor (%%) ...\n\n')
    fprintf(1,'         [%+.4f  %+.4f  %+.4f]\n',eps_print(1,1),eps_print(1,2),eps_print(1,3))
    fprintf(1,'   eps = [%+.4f  %+.4f  %+.4f]\n',eps_print(2,1),eps_print(2,2),eps_print(2,3))
    fprintf(1,'         [%+.4f  %+.4f  %+.4f]\n',eps_print(3,1),eps_print(3,2),eps_print(3,3))
    fprintf(1,'\n')
    fprintf(1,'and internal strain\n\n')
    fprintf(1,'   xi  = %.3f \n\n', par.xi)




  %%%%%%%%%%%%%%%%%%
  % diamond (reciprocal) lattice vectors (including strain)
    fprintf(1,'  compute (reciprocal) lattice vectors ... ')
    tic
    [par.a, par.b, par.tau] = diamond(par.a0, par.eps, par.xi);

  % create reciprocal lattice vectors
    [par.G_list, par.N_G] = reciprocal_lattice_vectors(par);
    toc


  %%%%%%%%%%%%%%%%
  % conduction band minimum
    fprintf(1,'  determine conduction band minimum ...... ')
    tic
    K_Delta      = find_minimum([0 0 0.85]', par); % search for k-point minimum (in units of 2pi/a0) in Si conduction band
    par.k0       = K_Delta(3) * 2*pi/par.a0;       % z-component of conduction band minimum k-vector 
    par.k0_vec   = [0 0 par.k0]';                  % store k-vector of conduction band minimum
    toc

  % k1 value
    par.k1       = 2*pi/par.a0*(1 - par.eps(3,3)) - par.k0;

  % compute effective mass  
    fprintf(1,'  compute effective mass ................. ')
    tic
    [mass_Delta] = compute_effective_mass(par.idx_CB_low, par.k0_vec, 1E-4 * 2*pi/par.a0, par); % finite difference approx of effective mass tensor at conduction band minimum
    par.mass_t   = 0.5*(mass_Delta(1,1)+mass_Delta(2,2)) * par.const.m0; % in SI units   
    par.mass_l   = mass_Delta(3,3) * par.const.m0; % in SI units
    toc

  % compute indirect band gap 
    fprintf(1,'  indirect band gap ...................... ')
    tic
    H = Hamiltonian(par.k0_vec, par);
    [E] = eig(H,'vector');
    Ec = E(par.idx_CB_low) * par.energy_scale;
    
    H = Hamiltonian([0 0 0]', par);
    [E] = eig(H,'vector');
    Ev = E(par.idx_VB_high) * par.energy_scale;
    
    par.Eg = Ec - Ev;
    toc
  
  % print results  
    fprintf(1,'\n')
    fprintf(1,'  Summary of conduction band parameter computation:\n')
    fprintf(1,'  ------------------------------------------------ \n')
    fprintf(1,'  conduction band minimum:      k0 = (%.4f, %.4f, %.4f) * 2*pi/a_0\n',K_Delta(1), K_Delta(2), K_Delta(3));
    fprintf(1,'  indirect band gap:            Eg = %.4f eV\n',par.Eg/par.units.eV);
    fprintf(1,'  effective mass tensor (m0):        [%+.4f  %+.4f  %+.4f],\n',mass_Delta(1,1),mass_Delta(1,2),mass_Delta(1,3));
    fprintf(1,'                                 m = [%+.4f  %+.4f  %+.4f],\n',mass_Delta(2,1),mass_Delta(2,2),mass_Delta(2,3));
    fprintf(1,'                                     [%+.4f  %+.4f  %+.4f],\n',mass_Delta(3,1),mass_Delta(3,2),mass_Delta(3,3));
    
  % compute band structure coefficients
    par.n_range = [-10:10];
    %[par.c, par.C2, par.C4, par.C4p] = compute_bandstructure_coefficients(par.n_range, par);
    [par.c, par.C2, par.C4] = compute_bandstructure_coefficients(par.n_range, par);


end

