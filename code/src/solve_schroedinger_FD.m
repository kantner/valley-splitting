function [S, E] = solve_schroedinger_FD(x, par)
% solve Schrödinger eigenvalue problem in longitudinal direction (growth direction)
% for prescribed epitaxial profile X = X_QW + x

  % grid
    N  = par.N;
    dz = par.dz;

  % number of eigenvalues
    neigs = par.neigs;

  % finite difference laplacian
    Lap = Laplacian(N);
    
  % kinetic energy operator
    T = - par.const.hbar*par.const.hbar/(2*par.mass_l) * 1/(dz^2) * Lap;

  % potential energy operator
    U_x = potential_modification(x, par);
    U   = par.U_QW + par.U_F + U_x;

  % find reference energy
    threshold = 0.99;    
    %idx       = find(par.QW_indicator > threshold);
    idx = find(par.QW_indicator > threshold * max(par.QW_indicator));
    %if isempty(idx)
    %  U_min = 0;
    %else
      U_min = min(U(idx));
    %end
    %U_min

  % Hamiltonian    
    H = (T + diag(sparse(U)) - U_min * speye(N))/par.energy_scale;
  
  % make Hermitian
    %H = 0.5*(H + H');

  % solve eigenvalue problem
    %opts.issym = 1;
    opts.disp  = 0;
    %opts.p     = 30;
     
    [S, E] = eigs(H, neigs, 'smallestabs', opts);


  % rescale energy to SI units  
    E = diag(E) * par.energy_scale + U_min;

  % DEBUG
    DEBUG = 0;
    if DEBUG == 1
      figure(243423);clf;hold all;
      plot(par.z/par.units.nm, U/par.units.meV, 'k-')
      for i = 1 : neigs
      %plot(par.z/par.units.nm, E(i)/par.units.meV + par.N * abs(S(:,i)).^2, 'b-')
      plot(par.z/par.units.nm, E(i)/par.units.meV + par.N * 10* U/par.units.eV .* abs(S(:,i)).^2, 'b-')
      end
      box on
      xlabel('z (nm)')
      ylabel('energy (meV)')
      title('DEBUG: Schrödinger solver')

      drawnow
      pause

    end


end

