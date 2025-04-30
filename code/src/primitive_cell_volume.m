function [Omega_p] = primitive_cell_volume(eps, par)

  % unstrained lattice vectors
    a = par.pp.a;

  % strained lattice vectors
    for i = 1 : 3
      a{i} = (eye(3) + eps) * a{i};
    end

  % primitive unit cell volume
    Omega_p = a{1}' * cross(a{2} , a{3});

end