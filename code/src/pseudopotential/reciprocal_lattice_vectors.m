function [par] = reciprocal_lattice_vectors(eps, par)

  % unstrained lattice vectors
    a = par.pp.a;

  % strained lattice vectors
    for i = 1 : 3
      a{i} = (eye(3) + eps) * a{i};
    end

  % primitive unit cell volume
    par.Omega_p = a{1}' * cross(a{2} , a{3});

  % atomic volume    
    par.Omega_a = 0.5 * par.Omega_p;

  % reciprocal lattice vectors
    par.pp.b{1} = 2*pi/par.Omega_p * cross(a{2} , a{3});
    par.pp.b{2} = 2*pi/par.Omega_p * cross(a{3} , a{1});
    par.pp.b{3} = 2*pi/par.Omega_p * cross(a{1} , a{2});


%% generate list of G-vectors   
   G_list   = [];
   k_max    = 0*[1,0,0]' * 2*pi/par.a0;
   G_max_sq = (2*par.const.m0)/par.const.hbar^2 * par.pp.E_cutoff;
   
 % set max integer for search sweep
   n_max = ceil(sqrt(G_max_sq)*(par.a0/(2*pi)));

   for nx = -n_max : n_max  
     G1 = nx * par.pp.b{1};
     for ny = -n_max : n_max
       G2 = G1 + ny * par.pp.b{2};
       for nz = -n_max : n_max
         G3 = G2 + nz * par.pp.b{3};

       % check cutoff criterion and append to list       
         if (G3 + k_max)'*(G3 + k_max) < G_max_sq
           G_list = [G_list, G3];
         end

       end
     end
   end

   par.pp.G_list = G_list;

   par.pp.N_G = size(par.pp.G_list,2);

end