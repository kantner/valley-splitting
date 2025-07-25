clear all
clc

 %% path
    addpath('colormaps/')
    addpath('src/')
    addpath('src/pseudopotential/')
    addpath('scripts/')
    addpath('scripts/silicon_lattice/')
    addpath('scripts/silicon_lattice/arrow3D_pub/')
    addpath('scripts/silicon_lattice/mArrow3')

 %% parameters
    physical_constants
    
  % constants  
    par.const.hbar = reducedPlanckConstant;
    par.const.e0   = elementaryCharge;
    par.const.m0   = electronMass;
    par.const.c0   = vacuumSpeedOfLight;
    par.const.aB   = bohrRadius;
    
  % units  
    par.units.eV   = elementaryCharge;
    par.units.meV  = 1E-3 * par.units.eV;
    par.units.ueV  = 1E-6 * par.units.eV;
    par.units.Ry   = 0.5*hartree;
    par.units.GPa  = 1E9;
    par.units.nm   = 1E-9;   
    par.units.Angstrom = 0.1 * par.units.nm;    
    par.units.percent  = 1E-2;

  % lattice constant
    par.a0_Si = 5.43 * par.units.Angstrom; % Fischetti % Laux (1996)
    par.a0_Ge = 5.65 * par.units.Angstrom; % Fischetti % Laux (1996)
    
    par.a0_SiGe = @(x) par.a0_Si + 0.200326 * par.units.Angstrom * x .* (1-x) + (par.a0_Ge-par.a0_Si)*x.^2; % Rieger % Vogl (1993)
    
    par.a0      = par.a0_Si;
    par.a0_subs = par.a0_SiGe(0.3);

  % add ML as a unit  
    par.units.ML = par.a0/4;
   
  % elasticity 
    par.C11 = 167.5 * par.units.GPa; % Rieger % Vogl (1993)
    par.C12 =  65.0 * par.units.GPa; % Rieger % Vogl (1993)
    par.C44 =  79.6 * par.units.GPa; % Fischetti % Laux (1996)
    
  % pseudopotential
    par.pp.nonlocal = 1; % 0 = off | 1 = on
    par.pp.SOI      = 0; % 0 = off | 1 = on (without SOI, the Ungersboeck pseudopotential model is identical to the model by Rieger & Vogl)

  % Fermi vector and energy
    par.pp.kF = 1.66 * 2*pi/par.a0; % Ungersboeck (2007)
    par.pp.EF = par.const.hbar^2 * par.pp.kF^2 /(2*par.const.m0);

  % local pseudopotential
    par.pp.V_0      = -2/3 * par.pp.EF;       % Rieger % Vogl (1993)
    par.pp.V_sqrt3  = -0.2241 * par.units.Ry; % Rieger % Vogl (1993)
    par.pp.V_sqrt8  = +0.0520 * par.units.Ry; % Rieger % Vogl (1993)
    par.pp.V_sqrt11 = +0.0724 * par.units.Ry; % Rieger % Vogl (1993)

  % compute cubic spline interpolation parameters 
    par.pp.BC_mode = 1; % different BCs for cubic spline interpolation. 1: V'(0)=0 | 2: V''(0) = 0
    [par.pp.w] = compute_cubic_spline_weights(par);

  % plot atomic pseudopotential
    %plot_atomic_pseudopotential(par)
    
  % nonlocal pseudopotential
    par.pp.A0 = 0.03 * par.units.Ry;       % s-well depth
    par.pp.R0 = 1.06 * par.units.Angstrom; % s-well radius

  % spin-orbit interaction
    par.pp.mu   = 0.00023 * par.units.Ry;
    par.pp.zeta = 7.5589  * 1/par.units.Angstrom;

  % set scale for internal computations  
    par.energy_scale = par.units.Ry;

  % cut off energy  
    par.pp.E_cutoff = 12 * par.units.Ry; % this gives 181 plane waves
  
  % diamond lattice vectors (relaxed)
    par.pp.a{1} = 0.50 * par.a0 * [0 1 1]';
    par.pp.a{2} = 0.50 * par.a0 * [1 0 1]';
    par.pp.a{3} = 0.50 * par.a0 * [1 1 0]';
    par.pp.tau  = 0.25 * par.a0 * [1 1 1]';

  % Kleinman's internal strain parameter
    par.pp.internal_strain = 0.53;

  % strain tensor  
    eps    = zeros(3,3);

    eps_QW = strain_quantum_well(0.3, par);
    %eps    = eps + eps_QW;

  % compute band structure parameters
    [par] = compute_conduction_band_parameters(eps, par);

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% heterostructure part
  % plot band structure
    %plot_bandstructure(eps, par)

  %%%%%%%%%%%%%%%%%%
  % quantum dot parameters
    par.omega_x = 3.0 * par.units.meV /par.const.hbar;
    par.omega_y = 3.0 * par.units.meV /par.const.hbar;
    par.l_x     = sqrt(par.const.hbar/(par.mass_t*par.omega_x));
    par.l_y     = sqrt(par.const.hbar/(par.mass_t*par.omega_y));

  % Si/Ge conduction band offset
    par.dEc     = 0.5 * par.units.eV;

  %%%%%%%%%%%%%%%%%  
  % quantum well parameters
    par.h_QW      = 75  * par.units.ML; % QW height
    par.X_barrier = 0.3;                % Ge content in barrier
    par.sigma_u   = 0.5 * par.units.nm; % width of upper interface
    par.sigma_l   = 0.5 * par.units.nm; % width of lower interface

  % vertical electric field (in V/m)
    par.F       = 5.0 * 1E6;

  %%%%%%%%%%%%%%%%
  % grid
    par.N  = 2^12;
    par.L  = 4*par.h_QW;
    par.dz = par.L/par.N;
    par.z  = [0:par.N-1]'*par.dz - par.L/2;

  % k-space grid  
    par.k = [0 : par.N/2-1, -par.N/2:-1]'*2*pi/par.L;

  %%%%%%%%%%%%%%%%
  % quantum well indicator function
    %par.QW_indicator = 0.5 * tanh((par.z+0.5*par.h_QW)/par.sigma_l) + 0.5*tanh((-par.z+0.5*par.h_QW)/par.sigma_u);  % QW at [-h/2 h/2]
    par.QW_indicator = 0.5 * tanh((par.z+par.h_QW)/par.sigma_l) + 0.5 * tanh((-par.z)/par.sigma_u);                 % QW at [-h 0]

  % nominal alloy profile
    par.X_QW = par.X_barrier * (1 - par.QW_indicator);
    
  % QW potential from nominal (smoothed) profile
    par.U_QW = par.dEc * par.X_QW;

  % electric field potential
    par.U_F  = -par.const.e0 * par.F * par.z;    
    
  %%%%%%%%%%%%%%%%
  % numerics

  % discretization method (kinetic energy operator)
    par.discretization = 1; % 1 = finite difference | 2 = spectral

  % number of eigenvalues
    par.neigs = 15; %  number of eigenvalues   

  % dummy value for artificially negative covariance functions
    par.Gamma_failure = 1E99;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% plot figures for paper 1

  % Fig. 1 (Brillouin zone, strained lattices)
    Fig1a_brillouin_zone_fcc
    Fig1cde_plot_lattice
  
  % Fig. 2 (band structure coefficients vs shear strain)
    Fig2_bandstructure_coeffs_vs_strain(par)

  % Fig. 3 
    Fig3_interface_width_dependency(par)
  
  % Fig. 4 
    Fig4_wiggle_well_line_plots(par)
  
  % Fig. 5 
    Fig5a_wiggle_well_wavenumber_and_amplitude_vs_shear_strain(par)
    Fig5bcde_wiggle_well_wavenumber_vs_amplitude(par)
  
  % Fig. 6
    Fig6_brillouin_zone_sectors
