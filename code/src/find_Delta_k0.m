function [k_Delta] = find_Delta_k0(k0_low, k0_up, tol_eV, par)



  k(:,1) = 2*pi/par.a0 * [0 0 k0_up]';
  k(:,2) = 2*pi/par.a0 * [0 0 k0_low]';

  k(:,3) = 0.5 *(k(:,1) + k(:,2));

  % initial data
    for i = 1 : 3
      H = Hamiltonian(k(:,i),par);
      E = eig(H,'vector');
      energy(i) = E(par.idx_CB_low);
    end

    maxiter = 100;

  
  continue_iteration = 1;

  while  continue_iteration == 1

  % pick lowest three points
    [~,idx]=sort(energy,'ascend');

  % quadratic fit    
    A = zeros(3,3);
    b = zeros(3,1);
    for i = 1 : 3
      for j = 1 : 3
        A(i,j) = (par.a0/(2*pi) * k(3,idx(i)) ).^(3-j);
      end
      b(i) = energy(idx(i));
    end

    coeff = A\b;

  % predict minimum
    k0_new = - 0.5*coeff(2)/coeff(1);

  % create new vector
    m = length(energy)+1;
    k(:,m) = 2*pi/par.a0 * [0 0 k0_new]';

  % compute energy
    H = Hamiltonian(k(:,m),par);
    E = eig(H,'vector');
    energy(m) = E(par.idx_CB_low);

  % check for termination
    if abs(energy(m) - energy(m-1)) < tol_eV
      continue_iteration = 0;
    end

    if m > maxiter
      continue_iteration = 0;
    end


  end


  k_Delta = k(:,end);

end


