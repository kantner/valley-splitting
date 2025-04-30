function [] = Fig7_confinement_dependency(par)
% sweep over interface width

  % curoff
    par.pp.E_cutoff = 6*par.units.Ry;

  % QW
    par.h_QW      = 75 * par.units.ML;
    par.X_barrier = 0.3;
 
  % grid
    par.N     = 2^13;
    par.L     = 5*par.h_QW;
    par.dz    = par.L/par.N;
    par.z     = [0:par.N-1]'*par.dz - par.L/2;

  % k-space grid  
    par.k = [0 : par.N/2-1, -par.N/2:-1]'*2*pi/par.L;

  % set strain
    [eps_QW] = strain_quantum_well(0.3, par);
    eps_xy   = 0.1 * par.units.percent;

    eps      = eps_QW;
    eps(1,2) = eps_xy;
    eps(2,1) = eps_xy;

  % quantum dot parameters
    par.omega_x = 3.0 * par.units.meV /par.const.hbar;
    par.omega_y = 3.0 * par.units.meV /par.const.hbar;
    par.l_x     = sqrt(par.const.hbar/(par.mass_t*par.omega_x));
    par.l_y     = sqrt(par.const.hbar/(par.mass_t*par.omega_y));
  
  % band structure parameters  
    [par] = compute_conduction_band_parameters(eps, par);

  % quantum well interface width
    par.sigma_u   = 0.5 * par.units.nm;
    par.sigma_l   = 0.5 * par.units.nm;
  
  %%%%%%%%%%%%%%%%
  % quantum well indicator
    z_mid = -0.5*par.h_QW;
    par.QW_indicator = 0.5 * tanh((par.z+0.5*par.h_QW)/par.sigma_l) + 0.5 * tanh((-par.z+0.5*par.h_QW)/par.sigma_u);
    par.QW_indicator = 0.5 * tanh((par.z+0.5*par.h_QW-z_mid)/par.sigma_l) + 0.5 * tanh((-par.z+0.5*par.h_QW+z_mid)/par.sigma_u);

  % nominal alloy profile    
    par.X_QW         = par.X_barrier * (1 - par.QW_indicator);
    
  % set electric field    
    par.F   = 0;
  
  % potential
    par.U_F = -par.const.e0 * par.F * par.z;     
 
  % wiggle well profile          
    x = zeros(par.N,1);

  % field strength
    dEc_range = logspace(-9,0,101)*500*par.units.meV;

  % sweep
    mean_E_VS = zeros(length(dEc_range),1);
    var_E_VS  = zeros(length(dEc_range),1);
    nu        = zeros(length(dEc_range),1);
    sigma     = zeros(length(dEc_range),1);
    Delta_det = zeros(length(dEc_range),1);

    for iEc = 1 : length(dEc_range)

    % set value  
      par.dEc = dEc_range(iEc);

    % QW potential from nominal (smoothed) profile  
      par.U_QW =1.17*par.units.eV + par.dEc * par.X_QW;      
    
    % compute valley splitting
      compute_derivatives = 0;
      [out] = compute_valley_splitting(x, eps, compute_derivatives, par);
      
      mean_E_VS(iEc) = out.M;
      var_E_VS(iEc)  = out.V;
      nu(iEc)        = out.nu;
      sigma(iEc)     = out.sigma;
      Delta_det(iEc) = out.Delta_det;

    % plot
      figure(21345432);clf;hold all;
        plot(par.z/par.units.nm, (par.U_QW + par.U_F)/par.units.meV, 'k-','LineWidth',2)
        for i = 1 : par.neigs
          plot(par.z/par.units.nm, out.E(i)/par.units.meV + par.N * abs(out.S(:,i)).^2,'-','Color',[1 1 1]*0.5)        
        end
        plot(par.z/par.units.nm, out.E(out.idx_gnd)/par.units.meV + par.N * abs(out.S(:,out.idx_gnd)).^2,'r-','LineWidth',2)        
        xlabel('space (nm)')
        ylabel('energy (meV)')
        box on
        drawnow

    end
    
  %%%%%%%%%%%%%%%
  % plot results
    figure(213454);clf;hold all;
      plot(dEc_range/par.units.meV , mean_E_VS/par.units.ueV, 'k-','LineWidth',2,'DisplayName','mean E_{VS}')
      plot(dEc_range/par.units.meV , sqrt(var_E_VS)/par.units.ueV, 'r-','LineWidth',2,'DisplayName','std E_{VS}')
      plot(dEc_range/par.units.meV , nu/par.units.ueV, 'b-','LineWidth',2,'DisplayName','\nu')
      plot(dEc_range/par.units.meV , sigma/par.units.ueV, 'm-','LineWidth',2,'DisplayName','\sigma')
      box on
      xlabel('\Delta E_c (meV)')
      ylabel('energy (ueV)')
      legend();



end
