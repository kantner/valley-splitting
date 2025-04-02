function [S, E] = solve_schroedinger_spectral(x, par)
% solve Schrödinger eigenvalue problem in longitudinal direction (growth direction)
% for prescribed epitaxial profile X = X_QW + x
% using spectral discretization

  % size
    N  = par.N;


  % number of eigenvalues
    neigs = par.neigs;

        
  % solve eigenvalue problem
    if par.neigs == N
      %opts.p     = min(20*par.neigs,par.N)      
    else
      opts.p     = 5*neigs%;min(10*par.neigs,par.N);     
    end
      opts.issym = 1;
      opts.disp  = 1;
      opts.maxit = 10000;

      [S, E] = eigs(@(psi) apply_H(psi, x, par), N, [], neigs, 'smallestreal',opts);




  % rescale energy to SI units  
    E = diag(E) * par.energy_scale;


  %%%%%%%%%%%%%%%%%%%%
  % DEBUG
  %{
        figure(2343545);clf;hold all;
      plot(par.z/par.units.nm, U/par.units.meV,'k-','DisplayName','U(z)','LineWidth',2)
      yline(U_min/par.units.meV,'r--','DisplayName','U_{min}','LineWidth',2)

  
  
      for i = 1 : neigs
        %yline(E(i)/par.units.meV,'b-')
        plot(par.z/par.units.nm, E(i)/par.units.meV + par.N * abs(S(:,i)).^2,'b-','DisplayName',['E_{',num2str(i),'}'],'LineWidth',2)
      end
      box on
      legend
      title('DEBUG: Schrödinger solver')
      
  drawnow
  %}

end


