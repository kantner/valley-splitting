function [] = evaluate_valley_splitting(par)

k_Delta = find_Delta_k0(0.7,0.9,1E-16,par);


k_vec = k_Delta;
H = Hamiltonian(k_vec, par);
[c, E] = eig(H);

c = c(:,par.idx_CB_low);
E = diag(E)*par.energy_scale;

% 2k0 theory factor
C0 = 0;
for i = 1 : par.N_G
  C0 = C0 + c(i) * c(par.idx_minus(i));
end
C0


% 2k0 theory factor
C0 = 0;
for i = 1 : par.N_G


  Gi = par.G_list(:,i);
  idx = find( sum(abs(par.G_list + Gi)) == 0);

  C0 = C0 + c(i) * c(idx);
end
C0


% 2k0 theory factor
C0 = 0;
for n = 1 : length(par.idx_QW)

  i = par.idx_QW(n,1);
  j = par.idx_QW(n,2);


  Gi = par.G_list(:,i);
  Gj = par.G_list(:,j);


  if abs(Gi(3) +Gj(3))*par.a0/(2*pi) < 1E-12
  C0 = C0 + c(i) * c(j);
  end
end
C0



% 2k0 theory factor
C0 = 0;
for n = 1 : length(par.idx_QW)

  i = par.idx_QW(n,1);
  j = par.idx_QW(n,2);


  Gi = par.G_list(:,i);
  Gj = par.G_list(:,j);


  if abs(Gi(3) +Gj(3))*par.a0/(2*pi) < 1E-12
  C0 = C0 + c(i) * c(j);
  end
end
C0


% 2k0 theory factor --- Gauss
disp('C0 ... Gauss')
C0 = 0;
a_QD = 1E-9;
for i = 1 : par.N_G
  Gi = par.G_list(:,i);
  for j = 1 : par.N_G
    Gj = par.G_list(:,j);


    % Gauss      
      delta_x = exp(-((Gi(1) + Gj(1))*a_QD)^2 );
      delta_y = exp(-((Gi(2) + Gj(2))*a_QD)^2 );
      delta_z = exp(-((Gi(3) + Gj(3))*a_QD)^2 );

    % sum  
      C0 = C0 + delta_x * delta_y * delta_z * c(i)*c(j);

  end
end
C0



% Long wavelength wiggle well --- Gauss
disp('C2 ... Gauss')
C2 = 0;
a_QD = 1E-9;
for i = 1 : par.N_G
  Gi = par.G_list(:,i);
  for j = 1 : par.N_G
    Gj = par.G_list(:,j);


    % Gauss      
      delta_x = exp(-((Gi(1) + Gj(1))*a_QD)^2 );
      delta_y = exp(-((Gi(2) + Gj(2))*a_QD)^2 );
      delta_z = exp(-((Gi(3) + Gj(3) - 4*pi/par.a0)*a_QD)^2 );

    % sum  
      new = delta_x * delta_y * delta_z * c(i)*c(j);
      %{
      if abs(new)>1E-8
        Gi*par.a0/(2*pi)
        Gj*par.a0/(2*pi)
        new
      end
      %}
      C2 = C2 + new;

  end
end
C2



% Long wavelength wiggle well --- Gauss
%
  Nq = 1001;
  q_range = linspace(0,2,Nq) * 2*pi/par.a0;
  Delta   = zeros(Nq,1);

  a_QD = 10E-9;
  k0   = 0.85 * 2*pi/par.a0;

  k0   = k_Delta(3);

  L    = 100E-9;
  
  for iq = 1 : Nq
    q = q_range(iq);
    for i = 1 : par.N_G
      Gi = par.G_list(:,i);
      for j = 1 : par.N_G
        Gj = par.G_list(:,j);
    
    
        % Gauss      
          delta_x = exp(-((Gi(1) + Gj(1))*a_QD)^2 );
          delta_y = exp(-((Gi(2) + Gj(2))*a_QD)^2 );
          
    
        % sum  
          Delta(iq) = Delta(iq) ...
                  + delta_x * delta_y * c(i)*c(j) * exp(0.5*1i*(Gi(3)+Gj(3)-2*k0)*L) * ( ...
                      exp(+1i*q*L/2) * sinc( (Gi(3)+Gj(3)-2*k0+q)*L/2 ) ...
                    + exp(-1i*q*L/2) * sinc( (Gi(3)+Gj(3)-2*k0-q)*L/2 ) ...
                    );
    
      end
    end
  end

  figure(123453);clf;hold all;
  plot(q_range * par.a0/(2*pi),abs(Delta),'k-')
  %plot(q_range * par.a0/(2*pi),real(Delta),'rx-')
  %plot(q_range * par.a0/(2*pi),imag(Delta),'bx-')
  set(gca,'YScale','log')
  box on
  title('cos(qz) contribution to valley splitting')
  xlabel('q (2\pi/a_0)')
  ylabel('weight')
  xline(2*k0 *par.a0/(2*pi),'r--')
  xline(2*(2*pi/par.a0-k0) *par.a0/(2*pi),'r--')


%}




end