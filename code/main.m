clear all
clc

format shorte


  %%%%%%%%%%%
  % add folders to path
    addpath('colormaps/')    
    addpath('scripts/')    
    addpath('scripts/silicon_lattice/')
    addpath('scripts/silicon_lattice/arrow3D_pub/')
    addpath('src/')

  %%%%%%%%%%%%%%%%%%%%
  % units and constants
    physical_constants;
    par.const.e0   = elementaryCharge;
    par.const.hbar = reducedPlanckConstant;
    par.const.m0   = electronMass;
    par.const.aB   = bohrRadius;

    par.units.nm   = 1E-9;
    par.units.eV   = elementaryCharge;
    par.units.meV  = 1E-3 * par.units.eV;
    par.units.ueV  = 1E-6 * par.units.eV;
    par.units.Ry   = rydberg;
    par.units.GPa  = 1E9;    

    par.energy_scale = par.units.eV;
  
  %%%%%%%%%%%%%%%%%%%%  
  % lattice constant
    par.a0_Si   = 0.543 * par.units.nm; % Si (Fischetti & Laux, 1996)
    par.a0_Ge   = 0.565 * par.units.nm; % Ge (Fischetti & Laux, 1996)
    par.a0_SiGe = @(x) (1-x) * par.a0_Si + x * par.a0_Ge - 1.88E-3*par.units.nm *x*(1-x); % R. A. Logan, J. M. Rowell, and F. A. Trumbore, Phys. Rev. 136, A1751 (1964)
    
  % Si lattice constant
    par.a0      = par.a0_Si;

  % atomic monolayer width  
    par.units.ML = par.a0/4;

  % atomic volume and volume of primitive cell
    par.Omega_a  = (0.5*par.a0)^3;  % atomic volume
    par.Omega_p  = 2 * par.Omega_a;  % primitive unit cell volume

  % Fermi wave number  
    par.KF       = (3*pi^2 * 4 / par.Omega_a)^(1/3); % cf. J. R. Chelikowsky and M. L. Cohen, Phys. Rev. B 10, 5095 (1974)    

  %%%%%%%%%%%%%%%%%%%%  
  % elastic constants for Si (Fischetti & Laux, 1996)
    par.C11 = 165.77 * par.units.GPa;
    par.C12 =  63.93 * par.units.GPa;
    par.C44 =  79.62 * par.units.GPa;

  %%%%%%%%%%%%%%%%%%%%  
  % strain tensor and internal strain parameter
    [par.eps_QW] = strain_quantum_well(0.3, par);
   
    par.eps = zeros(3,3);

  % biaxial strain (QW)    
    par.eps = par.eps_QW;
  
  % add shear strain
  %{
    par.eps(1,2) = 0.1 * 1E-2;
    par.eps(2,1) = par.eps(1,2);
  %}

  % internal ionic displacement parameter strain  
    par.xi  = 0.53;
  
  %%%%%%%%%%%%%%%%%%%%   
  % set special band indices
    par.idx_VB_high = 4; % highest valence band
    par.idx_CB_low  = 5; % lowest conduction band
  
  %%%%%%%%%%%%%%%%%%
  % pseudpotential parameters
  % cutoff energy (determines number of plane waves in pseudopotential calculation)  
    par.E_cutoff = 8.0 * par.units.Ry;
  
  % local form factors are encoded in the atomic potential (function V_atomic_Fischetti)
    par.V_atomic_model = 1; % 1 = Fischetti/Laux | 2 = Kim/Fischetti | 3 = Fischetti/Higman | 4 = Friedel  
  
  % non-local pseudopotential parameters (Fischetti & Laux, 1996)
    par.alpha0 = 0.55 * par.units.Ry; % s-well depth
    par.beta0  = 0.32;                % s-well energy dependence (unitless)
    par.R0     = 2.0 * par.const.aB;  % s-well radius

  %%%%%%%%%%%%%%%%%  
  % compute conduction band parameters
    [par] = compute_conduction_band_parameters(par.eps, par.xi, par);

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
    par.h_QW      = 75  * par.units.ML;
    par.X_barrier = 0.3;
    par.sigma_u   = 0.5 * par.units.nm;
    par.sigma_l   = 0.5 * par.units.nm;

  % electric field (in V/m)
    par.F       = 5.0 * 1E6;

  %%%%%%%%%%%%%%%%
  % grid
    par.N  = 2^11;
    par.L  = 5*par.h_QW;
    par.dz = par.L/par.N;
    par.z  = [0:par.N-1]'*par.dz - par.L/2;

  % k-space grid  
    par.k = [0 : par.N/2-1, -par.N/2:-1]'*2*pi/par.L;

  %%%%%%%%%%%%%%%%
  % quantum well indicator function
    par.QW_indicator = 0.5 * tanh((par.z+par.h_QW)/par.sigma_l) + 0.5 * tanh(-par.z/par.sigma_u);

  % nominal alloy profile    
    par.X_QW = par.X_barrier * (1 - par.QW_indicator);
    
  % QW potential from nominal (smoothed) profile  
    par.U_QW = par.dEc * par.X_QW;

  % electric field potential
    par.U_F  = -par.const.e0 * par.F * par.z;    
    
  %%%%%%%%%%%%%%%%
  % numerics

  % discretization method
    par.discretization = 1; % 1 = finite difference | 2 = spectral

  % number of eigenvalues
    par.neigs = 25; %  number of eigenvalues
    

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


  
     

