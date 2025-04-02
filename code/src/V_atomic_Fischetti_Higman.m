function [V] = V_atomic_Fischetti_Higman(G, par)
% input
%   G ..... scalar or list of scalars (modulus of lattice vector)
%   par ... struct containing system parameters
%
% output
%   V ..... atomic pseudpotential

% reference
%   Fischetti, M. V. / Higman, J. M. (1991)
%   Hess, Karl (Ed.)
%   Monte Carlo Device Simulation Theory and calculation of the deformation potential electron-phonon scattering rates in semiconductors


% check for negative values (V_atomic depends on non-negative abs(G) only)
  if sum(G<0) > 0
    error('G must not contain negative values.')
  end

% import some values
  aB = par.const.aB;

% convert to atomic units (length in Bohr radius)
  G = G * aB;

% pseudopotential parameters (table 2 in reference)
  a1 = -0.992;
  a2 =  0.791;
  a3 = -0.352;
  a4 = -0.18;
  %a5 = 5.0;
  %a6 = 0.3;

% atomic volume (in atomic units) 
  %Omega_a =  par.Omega_a / aB^3;

%%%%%%%%%%%%%%%%%%  
% effective potential
% basic curve
  V = a1./(G.^2) .* (cos(a2 * G) + a3) .* exp(a4 * G.^4);

% 
% add units (Rydberg) to arrive at SI base units
  V = V * par.units.Ry;  
%{
% tanh-cutoff function
  cutoff = 0.5* (1 + tanh( (a5 - G.^2)./a6) );

% multiply basic curve and cutoff 
  V = V0 .* cutoff;

%}
  

end