function [G_list, N_G] = reciprocal_lattice_vectors(par)

  %G_max = (2*pi)/par.a0 * sqrt(nmax);

  E_cutoff = par.E_cutoff;

  hbar = par.const.hbar;
  m0   = par.const.m0;
  a0   = par.a0;

  for i = 1 : 3
  b_rescaled{i} = a0/(2*pi) * par.b{i};
  end



  gmax_squared = 2*m0*E_cutoff * (a0/(2*pi*hbar))^2;

  nmax = ceil(sqrt(gmax_squared));

  g_list = [];
  for n1 = -nmax : nmax
    for n2 = -nmax : nmax
      for n3 = -nmax : nmax
        g = n1*b_rescaled {1} + n2*b_rescaled{2} + n3*b_rescaled{3};
        if  g'*g <= gmax_squared
        % append  
          g_list = [g_list, g];
        end
      end
    end
  end


  G_list = g_list*2*pi/a0;
  N_G = size(g_list,2);

end