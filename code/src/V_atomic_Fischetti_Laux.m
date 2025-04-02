function [V] = V_atomic_Fischetti_Laux(G, par)
% input
%   G ..... scalar or list of scalars (modulus of lattice vector)
%   par ... struct containing system parameters
%
% output
%   V ..... atomic pseudpotential

% reference
%   M. V. Fischetti; S. E. Laux, J. Appl. Phys. 80, 2234–2252 (1996)
%   DOI: 10.1063/1.363052


% check for negative values (V_atomic depends on non-negative abs(G) only)
  if sum(G<0) > 0
    error('G must not contain negative values.')
  end

% import some values
  aB = par.const.aB;

% convert to atomic units (length in Bohr radius)
  G = G * aB;

% pseudopotential parameters (table 2 in reference)
  a1 = 8.659814;
  a2 = 3.505847;
  a3 = 0.7218257;
  a4 = 3.654774;
  a5 = 5.0;
  a6 = 0.3;

% atomic volume (in atomic units) 
  Omega_a =  par.Omega_a / aB^3;

%%%%%%%%%%%%%%%%%%  
% effective potential
% basic curve
  V0 = 1/Omega_a * a1 .* (G.^2 - a2)./(exp(a3 * (G.^2 - a4)) + 1);

% add units (hartree) to arrive at SI base units
  V0 = V0 * 2*par.units.Ry;  

% tanh-cutoff function
  cutoff = 0.5* (1 + tanh( (a5 - G.^2)./a6) );

% multiply basic curve and cutoff 
  V = V0 .* cutoff;


  

end