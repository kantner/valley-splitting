function [S, E] = solve_schroedinger(x, par)
% solve Schrödinger eigenvalue problem in longitudinal direction (growth direction)
% for prescribed epitaxial profile X = X_QW + x
% 

% actually accuary should be higher when using spectral method, but I have
% observed worse results in adjoint equation then ... find out why ...
%[S, E] = solve_schroedinger_spectral(x, par);
  
  
  switch(par.discretization)
    case 1 % finite differences
      [S, E] = solve_schroedinger_FD(x, par);
    case 2 % spectral method
      [S, E] = solve_schroedinger_spectral(x, par);
    otherwise
      error('not implemented')
  end


end

