function [par] = compute_conduction_band_parameters(eps, par)

  % high symmetry points (relaxed)
    k_L     = [0.5 0.5 0.5]' * 2*pi/par.a0;
    k_Gamma = [0 0 0]'       * 2*pi/par.a0;
    k_X     = [0 0 1]'       * 2*pi/par.a0;
    k_K     = [3/4, 3/4, 0]' * 2*pi/par.a0;        
    k_U     = [1/4, 1/4, 1]' * 2*pi/par.a0;
    k_W     = [1/2, 0, 1]'   * 2*pi/par.a0;

 % band indices  
   switch(par.pp.SOI)     
     case 0
        par.pp.idx_SO = 2;
        par.pp.idx_LH = 3;
        par.pp.idx_HH = 4;
        par.pp.idx_CB = 5;
      case 1
        par.pp.idx_SO = 4;
        par.pp.idx_LH = 6;
        par.pp.idx_HH = 8;
        par.pp.idx_CB = 10;
   end
   
  %%%%%%% 
  % top valence band energy at Gamma point (as reference energy)
    k = k_Gamma;
    H = Hamiltonian(k, eps, par);
    [~,E] = eig(H);
    E = diag(E) * par.energy_scale;
    Ev_offset = E(par.pp.idx_HH);

  %%%%%%%  
  % SO splitting  
    Delta_SO = E(par.pp.idx_HH) - E(par.pp.idx_SO);

  %%%%%%%%%%%%%%%%%%%  
  % hole masses
    %[m_HH_PT] = compute_effective_mass_perturbation_theory(par.pp.idx_HH, k_Gamma, eps, par);

  % finite difference approx  
  %{
    dk = 1E-5 * 2*pi/par.a0;

    [m_HH] = compute_effective_mass_finite_difference(par.pp.idx_HH, k_Gamma, dk, eps, par);
    [m_LH] = compute_effective_mass_finite_difference(par.pp.idx_LH, k_Gamma, dk, eps, par);
    [m_SO] = compute_effective_mass_finite_difference(par.pp.idx_SO, k_Gamma, dk, eps, par);
  %}

  %%%%%%%  
  % band gap at L point (w.r.t to valence band Gamma point)
    k = k_L;
    H = Hamiltonian(k, eps, par);
    [~,E] = eig(H);
    E = diag(E) * par.energy_scale;
    Eg_L = E(par.pp.idx_CB) - Ev_offset;

  
  %%%%%%%  
  % band gap at Gamma point (w.r.t to valence band Gamma point)
    k = k_Gamma;
    H = Hamiltonian(k_Gamma, eps, par);
    [~,E] = eig(H);
    E = diag(E) * par.energy_scale;
    Eg_Gamma = E(par.pp.idx_CB) - Ev_offset;

    
  %%%%%%%  
  % electron effective mass at L
    dk    = 1E-5 * 2*pi/par.a0;
    [m_L] = compute_effective_mass_finite_difference(par.pp.idx_CB, k_L, dk, eps, par);

    %[m2] = compute_effective_mass_perturbation_theory(par.pp.idx_CB, k_L, eps, par);
    %m2/par.const.m0    
  %}

  %%%%%%%  
  % find Delta minimum
    [k_Delta] = find_k_Delta(eps, par);

  % export
    par.k_Delta = k_Delta;
    par.k0      = k_Delta(3);
    par.k1      = 2*pi/par.a0 * (1 - eps(3,3)) - par.k0;
  
  % visual check if minimum is reasonable    
    plot_visual_check = 0;

    if plot_visual_check == 1
      k_scale_range = k_Delta(3) * par.a0/(2*pi) + linspace(-1, 1, 51)*1E-2;
      
      Ec_range = zeros(1,length(k_scale_range));
      for ik = 1 : length(k_scale_range)
        k = [0 0 k_scale_range(ik)]' * 2*pi/par.a0;
        H = Hamiltonian(k, eps, par);
        [~,E] = eig(H);
        E = diag(E)*par.energy_scale;
        Ec_range(ik) = E(par.pp.idx_CB);
      end
  
      figure(2);clf;hold all;
        plot(k_scale_range, Ec_range/par.units.eV, 'ko-')
        xline(k_Delta(3) * par.a0/(2*pi),'r-')
        box on
        xlabel('k (2\pi/a_0)')
        ylabel('energy (eV)')
        title('conduction band minimum position')
    end    

  %%%%%%%
  % band gap at Delta point (w.r.t to valence band Gamma point)
    k = k_Delta;
    H = Hamiltonian(k, eps, par);
    [c,E] = eig(H);
    E = diag(E) * par.energy_scale;
    Eg_Delta = E(par.pp.idx_CB) - Ev_offset;    

    

  %%%%%%%  
  % effective mass at Delta    
    dk        = 1E-5 * 2*pi/par.a0;
    [m_Delta] = compute_effective_mass_finite_difference(par.pp.idx_CB, k_Delta, dk, eps, par);

    %[m2] = compute_effective_mass_perturbation_theory(par.pp.idx_CB, k_Delta, eps, par);
    %m2/par.const.m0

    par.mass_l = m_Delta(3,3);
    par.mass_t = 0.5 * (m_Delta(1,1) + m_Delta(2,2) );
    
  %%%%%%%%%%%%%  
  % Bloch factor plane wave expansion coefficients and bandstructure coefficients C^(2), C^(4) for conduction band
    %par.c = c(:,par.pp.idx_CB);
    par.n_range = [-15:15]';
    [par.c, par.C2, par.C4] = compute_bandstructure_coefficients(par.n_range, eps, par);


  %%%%%%%%%%%%%%%%%%%   
  % cell volume   
    [par.Omega_p] = primitive_cell_volume(eps, par);

  % atomic volume  
    par.Omega_a   = 0.5 * par.Omega_p;


    
  %%%%%%%%%%%%%%    
  % conduction band deformation potential at Delta
    deps = 1E-8;
    [D] = compute_deformation_potential(par.pp.idx_CB, k_Delta, eps, deps, par);


  % valence band deformation potential at Gamma
    deps = 1E-8;
    [D] = compute_deformation_potential(par.pp.idx_HH, k_Gamma, eps, deps, par);

    

    fprintf('\n')
    fprintf('summary of band structure parameters ...\n')
    fprintf('  energy gap at L      Eg = %.4f eV\n',Eg_L/par.units.eV)
    fprintf('  energy gap at Delta  Eg = %.4f eV\n',Eg_Delta/par.units.eV)
    fprintf('  energy gap at Gamma  Eg = %.4f eV\n',Eg_Gamma/par.units.eV)
    fprintf('  Delta point:         k  = [%.4f, %.4f, %.4f] * 2 pi/a0\n',k_Delta(1)*par.a0/(2*pi),k_Delta(2)*par.a0/(2*pi), k_Delta(3)*par.a0/(2*pi));
    fprintf('  electron mass (Delta):\n')
    m = m_Delta/par.const.m0;
    fprintf('                            [%.3f,  %.3f,  %.3f]\n',m(1,1),m(1,2),m(1,3))
    fprintf('                       m  = [%.3f,  %.3f,  %.3f]\n',m(2,1),m(2,2),m(2,3))
    fprintf('                            [%.3f,  %.3f,  %.3f]\n',m(3,1),m(3,2),m(3,3))
    

end