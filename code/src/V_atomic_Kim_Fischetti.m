function [V] = V_atomic_Kim_Fischetti(G, par)
% input
%   G ... scalar or list of scalars (modulus of lattice vector)

% only non-negative values
  if sum(G<0) > 0
    error('G must not contain negative values.')
  end

% import
  a0  = par.a0;
  %V_0 = par.V_0;
  V_0  = -1.113*par.units.Ry;
  EF = -3/2 * V_0;
  kF = sqrt(2*par.const.m0*EF)/par.const.hbar;

% rescale
  g = G*a0/(2*pi);

% compute spline parametes  
  [c] = cubic_spline_parameters(par);
  par.c_spline = c;

% interval end points
  g_interval =  [0, 2*pi/a0*sqrt(3), 2*pi/a0*sqrt(8), 2*pi/a0*sqrt(11), 3*kF] * a0/(2*pi);

  idx1 = g_interval(1)<=g & g<g_interval(2);
  idx2 = g_interval(2)<=g & g<g_interval(3);
  idx3 = g_interval(3)<=g & g<g_interval(4);
  idx4 = g_interval(4)<=g & g<=g_interval(5);

  V = zeros(size(G));

  V(idx1) = par.c_spline(1) * g(idx1).^3 + par.c_spline(2) * g(idx1).^2 + par.c_spline(3) * g(idx1) + par.c_spline(4);
  V(idx2) = par.c_spline(5) * g(idx2).^3 + par.c_spline(6) * g(idx2).^2 + par.c_spline(7) * g(idx2) + par.c_spline(8);
  V(idx3) = par.c_spline(9) * g(idx3).^3 + par.c_spline(10) * g(idx3).^2 + par.c_spline(11) * g(idx3) + par.c_spline(12);
  V(idx4) = par.c_spline(13) * g(idx4).^3 + par.c_spline(14) * g(idx4).^2 + par.c_spline(15) * g(idx4) + par.c_spline(16);
  

% multiply with tanh-cutoff
  aB = par.const.aB;
  a5 = 5.0 /aB^2;
  a6 = 0.3 /aB^2;

  cutoff = 0.5* (1 + tanh( (a5 - G.^2)./a6) );

  V = V .* cutoff;

  

end