function [K_min] = find_minimum(K_vec0, par)
% find minimum of the conduction band near K_vec0


%{

% map 
  Nx = 55;
  Nz = 54;
  E = zeros(Nx,Nz);
  Kx = linspace(0,1,Nx);
  Kz = linspace(0,1,Nz);
  for ix = 1 : Nx
    ix
    for iz = 1 : Nz
      K_vec = [Kx(ix) 0 Kz(iz)]';
      E(ix,iz) =  conduction_band_energy(K_vec, par);       
    end
  end



figure(23435);clf; hold all;
  surf(Kx,Kz,E')

%}
  %return


% bounds
  K_ub = [+1 +1 +1];
  K_lb = [-1 -1 -1];


% opts
  opts = optimoptions(@fmincon, 'Display','notify-detailed','OptimalityTolerance',1E-16,'StepTolerance',1E-16,'ConstraintTolerance',1E-16,'FiniteDifferenceType','central');

  K_min = fmincon(@(K_vec) conduction_band_energy(K_vec, par), K_vec0,[],[],[],[],K_lb,K_ub,[],opts);

  

% fminunc for comparison
%{
  opts = optimoptions(@fminunc, 'Display','iter-detailed','OptimalityTolerance',1E-16,'StepTolerance',1E-16,'FiniteDifferenceType','central');
  K_min = fminunc(@(K_vec) conduction_band_energy(K_vec, par), K_vec0, opts)
%}
end