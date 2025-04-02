function [V,V0,cutoff] = V_atomic_Friedel(G, par)
% input
%   G ... scalar or list of scalars (modulus of lattice vector)

% only non-negative values
  if sum(G<0) > 0
    error('G must not contain negative values.')
  end

% import
  a0  = par.a0;

  aB = par.const.aB;

% Friedel parameters
  a1 = 106.0686;
  a2 = 2.2278;
  a3 = 0.6060;
  a4 = -1.9720;
  a5 = 5.0;
  a6 = 0.3;


% rescale
  %g = G*a0/(2*pi);

% convert to atomic units
  G = G*aB;

% basic curve
  V0 = a1 .* (G.^2 - a2)./(exp(a3 * (G.^2 - a4)) + 1);


% tanh-cutoff
  cutoff = 0.5* (1 + tanh( (a5 - G.^2)./a6) );

% multiply  
  V = V0 .* cutoff;

% add units (times hartree) to arrive at SI base units,
  V = V * 2*par.units.Ry;

% normalize by cell volume
  V = V * par.const.aB^3/(par.Omega_a);
  

end