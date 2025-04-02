function [] = plot_bandstructure(par)
  
  % import
    a0  = par.a0;
    N_G = par.N_G;

  %%%%%%%
  % high symmetry points
    k_Gamma   = [0, 0, 0]'       * 2*pi/a0;
    %k_X       = [0, 1, 0]'       * 2*pi/a0;
    k_Xx       = [1, 0, 0]'       * 2*pi/a0;    
    k_Xy       = [0, 1, 0]'       * 2*pi/a0;
    k_Xz       = [0, 0, 1]'       * 2*pi/a0;
    k_L       = [1/2, 1/2, 1/2]' * 2*pi/a0;
    k_K       = [3/4, 3/4, 0]'   * 2*pi/a0;    
    %k_U       = [1/4, 1, 1/4]'   * 2*pi/a0;
    k_U       = [1/4, 1/4, 1]'   * 2*pi/a0;
    k_W       = [1/2, 0, 1]'     * 2*pi/a0;


  % strained X points
    k_Xx       = 0.5*(par.b{2} + par.b{3});
    k_Xy       = 0.5*(par.b{1} + par.b{3});
    k_Xz       = 0.5*(par.b{1} + par.b{2});

    

%{
    k_Gamma = [0 0 0]';
    %k_X     = [0 0 1]';
    k_X     = [0 1 0]';
    k_L     = 0.5*[1 1 1]';
    
%}

  %%%%
  % create figure
    figure(2435); clf; hold all;
      box on
      xlabel('wave vector')
      ylabel('energy (eV)')

      ylim([-12.5, 6.5])
    
      title('silicon')


    N_pts = 101;
    s = linspace(0,1,N_pts);

  %%%%%%%%%%%%%%%%%%%%%
  % routes
    K_init{1}   = k_L;
    K_final{1}  = k_Gamma;

    K_init{2}   = k_Gamma;
    K_final{2}  = k_Xz;

    K_init{3}   = k_Xz;
    K_final{3}  = k_U;

    K_init{4}   = k_K;
    K_final{4}  = k_Gamma;

    K_init{5}   = k_Gamma;
    K_final{5}  = -k_Xz;

    for i = 1 : 5
      K_length(i) = norm(K_final{i}-K_init{i})*par.a0/(2*pi);
    end

    




  % compute valence band at Gamma for reference
    H = Hamiltonian(k_Gamma, par);
    E = eig(H,'vector');
    Ev = E(par.idx_VB_high) * par.energy_scale;

  %%%%%%
  % compute band diagram
    for iK = 1 : length(K_init)

    % assemble  
      energy = zeros(N_G, N_pts);
      for is = 1 : N_pts
        k_vec = (1-s(is)) * K_init{iK} + s(is) * K_final{iK};
  
        H = Hamiltonian(k_vec, par);
        E = eig(H,'vector');
        
        energy(:,is) = E * par.energy_scale;
      end
  
      for j = 1 : N_G
        switch(j)
          case par.idx_CB_low
            color = [1 0 0];
          case par.idx_VB_high 
            color = [0 0 1];  
          otherwise
            color = [0 0 0];
        end
      plot(sum(K_length(1:iK-1)) + K_length(iK)*s,(energy(j,:)-Ev)/par.units.eV,'-','Color',color,'LineWidth',2)
      end

    
      xline(sum(K_length(1:iK-1)),'k--')
    end


    set(gca,'XTick',[0,cumsum(K_length)])
    set(gca,'XTickLabel',{'L','\Gamma','X(z)','U;K','\Gamma','X(y)'})


    drawnow




  %%%%%%  
  % compute band gaps and masses
    dk = 1E-3*2*pi/par.a0; % finite difference

  % Delta
    %k_Delta = find_Delta_k0(0.0,1.0,1E-16,par);
    K_Delta = find_minimum([0 0 0.85]', par);
    k_Delta = 2*pi/par.a0 * K_Delta;

    k0 = k_Delta(3) * par.a0/(2*pi)


    k_vec = k_Delta;
    H = Hamiltonian(k_vec,par);
    [c,E] = eig(H);
    E = diag(E)*par.energy_scale;
    Eg_Delta = (E(par.idx_CB_low) - Ev)/par.units.eV
    yline(Eg_Delta,'r--')

    mass_Delta = compute_effective_mass(par.idx_CB_low, k_vec, dk, par)

%{
    c = c(:,par.idx_CB_low);

    g_list = par.G_list*a0/(2*pi);
    for i = 1 : par.N_G
    fprintf(1,'%d\t%d\t%d\t%.4e\n',round(g_list(1,i)),round(g_list(2,i)),round(g_list(3,i)),c(i) )
    end

    C0 = 0;
    for i = 1 : par.N_G

      Gi  = par.G_list(:,i);
      idx = find(sum(abs(par.G_list+Gi))==0);
      C0 = C0 + c(i) * c(idx);

    end
    C0
%}

  % L
    k_vec = k_L;
    H = Hamiltonian(k_vec,par);
    E = eig(H,'vector')*par.energy_scale;
    Eg_L = (E(par.idx_CB_low) - Ev)/par.units.eV

    mass_L = compute_effective_mass(par.idx_CB_low, k_vec, dk, par)

  % X(z)
    k_vec = k_Xz;
    H = Hamiltonian(k_vec,par);
    E = eig(H,'vector')*par.energy_scale;
    Eg_X = (E(par.idx_CB_low) - Ev)/par.units.eV

    mass_X = compute_effective_mass(par.idx_CB_low, k_vec, dk, par)


  % Gamma  
    k_vec = k_Gamma;
    H = Hamiltonian(k_vec,par);
    E = eig(H,'vector')*par.energy_scale;
    Eg_Gamma = (E(par.idx_CB_low) - Ev)/par.units.eV

    mass_Gamma = compute_effective_mass(par.idx_CB_low, k_vec, dk, par)
    mass_Gamma_holes = compute_effective_mass(par.idx_VB_high, k_vec, dk, par)




  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % create another figure
    figure(2436); clf; hold all;
  
    
      title('silicon')


    N_pts = 101;
    s = linspace(0,1,N_pts);

  %%%%%%%%%%%%%%%%%%%%%
  % routes
    clear('K_init','K_final')
    K_init{1}   =  k_Gamma;
    K_final{1}  =  k_Xx;

    K_init{2}   =  k_Gamma;
    K_final{2}  =  k_Xy;

    K_init{3}   =  k_Gamma;
    K_final{3}  =  k_Xz;



  %%%%%%%%%%
  % compute valence band at Gamma for reference
    H = Hamiltonian(k_Gamma, par);
    E = eig(H,'vector');
    Ev = E(par.idx_VB_high) * par.energy_scale;

  %%%%%%
  % compute band diagram
    for iK = 1 : length(K_init)
    subplot(1,3,iK); hold all;
      box on
      xlabel('wave vector')
      ylabel('energy (eV)')
      ylim([-12.5, 6.5])

    % assemble  
      energy = zeros(N_G, N_pts);
      for is = 1 : N_pts
        k_vec = (1-s(is)) * K_init{iK} + s(is) * K_final{iK};
  
        H = Hamiltonian(k_vec, par);
        E = eig(H,'vector');
        
        energy(:,is) = E * par.energy_scale;
      end
  
      for j = 1 : N_G
        switch(j)
          case par.idx_CB_low
            color = [1 0 0];
          case par.idx_VB_high 
            color = [0 0 1];  
          otherwise
            color = [0 0 0];
        end
      plot(s,(energy(j,:)-Ev)/par.units.eV,'-','Color',color,'LineWidth',2)
      end

    
      
    end





    drawnow

    
end