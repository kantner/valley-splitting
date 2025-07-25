function [] = Fig5bcde_wiggle_well_wavenumber_vs_amplitude(par)
% make plot like Woods (2024), Fig 3 (a)
%clear('map')

  % exportgraphics
    save_plot = 1; % 0 = off | 1 = on

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % strain tensor for QW
    [eps_QW] = strain_quantum_well(0.3, par);

  % do not compute derivates of valley-splitting functions
    compute_derivatives = 0;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Ge conc vs. q
    q_range = linspace(0,0.5,151) * 2*pi/par.a0;
    x_range = linspace(0,0.25,161);

  % set QW strain + fixed shear strain 
    %eps_xy = 0.1 * par.units.percent;
    eps_xy = 0.025 * par.units.percent;

    eps      = eps_QW;
    eps(1,2) = eps_xy;
    eps(2,1) = eps_xy;

  % compute conduction band parameters
    [par] = compute_conduction_band_parameters(eps, par);  

  %%%%%%%%%%%%%%%%%  
  % quantum well parameters
    par.h_QW      = 75  * par.units.ML;
    par.X_barrier = 0.3;
    par.sigma_u   = 0.5 * par.units.nm;
    par.sigma_l   = 0.5 * par.units.nm;

  % electric field (in V/m)
    par.F       = 5.0 * 1E6;

  %%%%%%%%%%%%%%%%
  % grid
    par.N  = 2^13;
    par.L  = 4*par.h_QW;
    par.dz = par.L/par.N;
    par.z  = [0:par.N-1]'*par.dz - 2*par.h_QW;

  % k-space grid  
    par.k = [0 : par.N/2-1, -par.N/2:-1]'*2*pi/par.L;
   %%
  %%%%%%%%%%%%%%%%
  % quantum well indicator
    %par.QW_indicator = 0.5 * tanh((par.z+0.5*par.h_QW)/par.sigma_l) + 0.5*tanh((-par.z+0.5*par.h_QW)/par.sigma_u); % QW at [-h/2 h/2]
    par.QW_indicator = 0.5 * tanh((par.z+par.h_QW)/par.sigma_l) + 0.5 * tanh((-par.z)/par.sigma_u); % QW at [-h 0]

  % nominal alloy profile    
    par.X_QW         = par.X_barrier * (1 - par.QW_indicator);
    
  % QW potential from nominal (smoothed) profile  
    par.U_QW = par.dEc * par.X_QW;

  % electric field potential
    par.U_F  = -par.const.e0 * par.F * par.z;    


  % allocate memory
    map.mean_E_VS = zeros(length(q_range),length(x_range));
    map.std_E_VS  = zeros(length(q_range),length(x_range));
    map.nu        = zeros(length(q_range),length(x_range));
    map.sigma     = zeros(length(q_range),length(x_range));
  
  % compute map
    fprintf(1,'compute map ... \n')
    tic
    for iq = 1 : length(q_range)

    % print progress  
      fprintf(1,'%.3f %%\n',iq/length(q_range)*100);

      for ix = 1 : length(x_range)
      
      % wiggle well profile          
        x_ww = x_wiggle_well(x_range(ix), q_range(iq), 0, par);

      % compute valley splitting
        [out] = compute_valley_splitting(x_ww, eps, compute_derivatives, par);
    
        mean_E_VS = out.M;
        var_E_VS  = out.V;
        nu        = out.nu;
        sigma     = out.sigma;

        map.mean_E_VS(iq, ix) = mean_E_VS;
        map.std_E_VS(iq, ix)  = sqrt(var_E_VS);
        map.nu(iq, ix)        = nu;
        map.sigma(iq, ix)     = sigma;

      end
    end
    toc


%% %%%%%%%%%%%%%%%
  % overview plot
    fig_obj = figure(90001); clf; hold all;
    sgtitle(['shear strain \epsilon_{x,y} = ',num2str(eps_xy * 100),' %'])

    nrows = 2;
    ncols = 3;

    cmap = divergingColormap1(256);
    xaxis_min = 0;
    xaxis_max = 0.5;
    yaxis_min = 0;
    yaxis_max = 0.2;


    subplot(nrows,ncols,1); hold all;
      surf(q_range*par.a0/(2*pi), x_range, (map.mean_E_VS/par.units.meV)')
      xlabel('q (2\pi/a_0)')
      ylabel('x (Ge amplitude)')
      zlabel('mean E_{VS} (meV)')      
      title(['mean(E_{VS})'])
      colorbar
      box on
      shading interp
      xlim([xaxis_min xaxis_max])
      ylim([yaxis_min yaxis_max])

    subplot(nrows,ncols,2); hold all;
      surf(q_range*par.a0/(2*pi), x_range, (map.std_E_VS/par.units.meV)')
      xlabel('q (2\pi/a_0)')
      ylabel('x (Ge amplitude)')
      zlabel('std E_{VS} (meV)')
      title(['std(E_{VS})'])
      colorbar
      box on
      shading interp
      colormap(cmap)
      xlim([xaxis_min xaxis_max])
      ylim([yaxis_min yaxis_max])

    subplot(nrows,ncols,3); hold all;
      surf(q_range*par.a0/(2*pi), x_range, (map.std_E_VS./map.mean_E_VS)')
      xlabel('q (2\pi/a_0)')
      ylabel('x (Ge amplitude)')
      zlabel('mean E_{VS} (meV)')      
      title(['ratio std(E_{VS})/mean(E_{VS})'])
      colorbar
      box on
      shading interp
      colormap(cmap)
      xlim([xaxis_min xaxis_max])
      ylim([yaxis_min yaxis_max])
  
    subplot(nrows,ncols,4); hold all;
      surf(q_range*par.a0/(2*pi), x_range, (map.nu/par.units.meV)')
      xlabel('q (2\pi/a_0)')
      ylabel('x (Ge amplitude)')
      zlabel('mean E_{VS} (meV)')
      title(['Rice parameter \nu'])
      colorbar    
      box on
      shading interp
      colormap(cmap)
      xlim([xaxis_min xaxis_max])
      ylim([yaxis_min yaxis_max])

    subplot(nrows,ncols,5); hold all;
      surf(q_range*par.a0/(2*pi), x_range, (map.sigma/par.units.meV)')
      xlabel('q (2\pi/a_0)')
      ylabel('x (Ge amplitude)')
      zlabel('mean E_{VS} (meV)')
      title(['Rice parameter \sigma'])
      colorbar   
      box on
      shading interp
      colormap(cmap)
      xlim([xaxis_min xaxis_max])
      ylim([yaxis_min yaxis_max])

    subplot(nrows,ncols,6); hold all;
      surf(q_range*par.a0/(2*pi), x_range, (map.nu./(0.3507 * 2*map.sigma))')
      xlabel('q (2\pi/a_0)')
      ylabel('x (Ge amplitude)')
      zlabel('\nu/(2\sigma)')
      title(['ratio of Rice parameters \nu/(2\sigma)'])
      colorbar   
      box on
      shading interp
      colormap(cmap)
      xlim([xaxis_min xaxis_max])
      ylim([yaxis_min yaxis_max])

  %%%%%%%%%%%%%
  % contour: deterministic enhancement separatrix
    contour_levels = [0.0 1.0]; % we need to add a second contour level, because somehow it does not do it for a single contour level alone (add 0 as dummy)
    M = contour(q_range*par.a0/(2*pi), x_range, (map.nu./(0.3507 * 2*map.sigma))',contour_levels);
     
  % total array size
    n = size(M,2);

    counter = 1;

  % plot all contour lines
    while counter <= n

      val   = M(1,counter);
      n_seg = M(2,counter); 
      x = M(1,counter+1: counter+n_seg);
      y = M(2,counter+1: counter+n_seg);
  
      if val == 1
        for i = 1 : nrows*ncols
          subplot(nrows,ncols,i); hold all;  
          z_dummy = 100*ones(size(y));
          plot3(x,y,z_dummy,'k-','LineWidth',2)
        end
      end

      counter = counter+n_seg+1;
    end


  %%%%%%%%%%%%%
  %{
  % contour: mean E_VS levels
    contour_levels = [0.5 1 2 3 4 5]*par.units.meV; % we need to add a second contour level, because somehow it does not do it for a single contour level alone (add 0 as dummy)
    M = contour(q_range*par.a0/(2*pi), x_range, map.mean_E_VS',contour_levels);
     
  % total array size
    n = size(M,2);

    counter = 1;

  % plot all contour lines
    while counter <= n

      val   = M(1,counter);
      n_seg = M(2,counter); 
      x = M(1,counter+1: counter+n_seg);
      y = M(2,counter+1: counter+n_seg);
  
      %if val == 1
        for i = 1 : nrows*ncols
          subplot(nrows,ncols,i); hold all;  
          z_dummy = 100*ones(size(y));
          plot3(x,y,z_dummy,'w-','LineWidth',2)
        end
      %end

      counter = counter+n_seg+1;
    end
  %}


  %%%%%%%%%%%
  % add lines        
    for i = 1 : nrows*ncols
      subplot(nrows,ncols,i); hold all;  
      xline(2*par.k0 * par.a0/(2*pi),'w--','LineWidth',2)
      xline(2*par.k1 * par.a0/(2*pi),'w--','LineWidth',2)
      xline(  par.k1 * par.a0/(2*pi),'m--','LineWidth',2)
    end
        
  % save
    drawnow
    if save_plot == 1
      height = 800;
      width  = 1600;
      set(fig_obj,'Position',[100 100 width height]);
      exportgraphics(fig_obj,'Fig5b_map_WW_E_VS_q_vs_x_at_eps_xy_fixed.pdf','ContentType','vector')
    end




%% %%%%%%%%%%%%%%%
  % plot for publication
    fig_obj = figure(90004); clf; hold all;
    sgtitle(['shear strain \epsilon_{x,y} = ',num2str(eps_xy * 100),' %'])

    nrows = 1;
    ncols = 1;

      %subplot(nrows,ncols,1); hold all;
      surf(q_range*par.a0/(2*pi), x_range, (map.mean_E_VS/par.units.ueV)')
      xlabel('q (2\pi/a_0)')
      ylabel('x (Ge amplitude)')
      zlabel('mean E_{VS} (meV)')      
      title(['mean(E_{VS})'])
      colorbar
      box on
      shading interp      
      xlim([0 0.5])
      ylim([0 0.25])
      cmap2 = divergingColormap1(16);
      %cmap2 = viridis(20);
      colormap(cmap2)
      clim([0,400])
      



  %%%%%%%%%%%%%
  % contour: deterministic enhancement separatrix
    contour_levels = [0.0 1.0]; % we need to add a second contour level, because somehow it does not do it for a single contour level alone (add 0 as dummy)
    M = contour(q_range*par.a0/(2*pi), x_range, (map.nu./(0.3507 * 2*map.sigma))',contour_levels);
     
  % total array size
    n = size(M,2);

    counter = 1;

  % plot all contour lines
    while counter <= n

      val   = M(1,counter);
      n_seg = M(2,counter); 
      x = M(1,counter+1: counter+n_seg);
      y = M(2,counter+1: counter+n_seg);
  
      if val == 1
        for i = 1 : nrows*ncols
          %subplot(nrows,ncols,i); hold all;  
          z_dummy = 5000*ones(size(y));
          plot3(x,y,z_dummy,'k-','LineWidth',2)
        end
      end

      counter = counter+n_seg+1;
    end


  %%%%%%%%%%%%%
  %{
  % contour: mean E_VS levels
    contour_levels = [0.5 1 2 3 4 5]*par.units.meV; % we need to add a second contour level, because somehow it does not do it for a single contour level alone (add 0 as dummy)
    M = contour(q_range*par.a0/(2*pi), x_range, map.mean_E_VS',contour_levels);
     
  % total array size
    n = size(M,2);

    counter = 1;

  % plot all contour lines
    while counter <= n

      val   = M(1,counter);
      n_seg = M(2,counter); 
      x = M(1,counter+1: counter+n_seg);
      y = M(2,counter+1: counter+n_seg);
  
      %if val == 1
        for i = 1 : nrows*ncols
          subplot(nrows,ncols,i); hold all;  
          z_dummy = 100*ones(size(y));
          plot3(x,y,z_dummy,'w-','LineWidth',2)
        end
      %end

      counter = counter+n_seg+1;
    end
  %}

  %%%%%%%%%%%
  % add lines    
    for i = 1 : nrows*ncols
      %subplot(nrows,ncols,i); hold all;  
      xline(2*par.k0 * par.a0/(2*pi),'w--','LineWidth',2)
      xline(2*par.k1 * par.a0/(2*pi),'w--','LineWidth',2)
      xline(  par.k1 * par.a0/(2*pi),'m--','LineWidth',2)
    end
        
  % save
    drawnow
    if save_plot == 1
      height = 800;
      width  = 1600;
      set(fig_obj,'Position',[100 100 width height]);
      exportgraphics(fig_obj,'Fig5b_map_WW_E_VS_q_vs_x_at_eps_xy_fixed_mean.pdf','ContentType','vector')
    end




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % make plots of wave function at different points
    fig_obj = figure(90003); clf; hold all;
    sgtitle(['shear strain \epsilon_{x,y} = ',num2str(eps_xy * 100),' %'])
  

    nrows = 5;
    ncols = 5;

    plot_counter = 0;

    for iq = 1 : 5
  
      switch(iq)
        case 1 % uniform Ge
          x = 0.05;
          q = 0;
        case 2 % "Ge spike" (harmonic/second order effect)
          x = 0.20;
          q = par.k1;
        case 3 % LP-WW
          x = 0.15;
          q = 2*par.k1;
        case 4 % harmonic/second order effect
          x = 0.05;
          q = par.k0;
        case 5 % SP-WW
          x = 0.05;
          q = 2*par.k0;
          
      end

      %{
      figure(90002); hold all;
        subplot(1,1,1); hold all;
        plot3(q*par.a0/(2*pi),x,100,'d','MarkerSize',10,'MarkerFaceColor','k')
      %}
        
      
      
    
      x_ww = x_wiggle_well(x, q, 0, par);   
  
      [out] = compute_valley_splitting(x_ww, eps, compute_derivatives, par);
      
      mean_E_VS   = out.M;
      std_E_VS    = sqrt(out.V);
      nu          = out.nu;
      sigma       = out.sigma; % this is the sigma for the rician (not for normal!)
      Gamma       = 0.5 * sigma^2;
      Delta_det   = out.Delta_det;
      Delta_det_n = out.Delta_det_n;
        
      plot_counter = plot_counter + 1;

      figure(90003); hold all;

      %  potential and wave function
      subplot(nrows,ncols, 0*ncols + iq); hold all;
        plot(par.z/par.units.nm, (par.U_F + par.U_QW)/par.units.meV, 'b-', 'LineWidth',1)
        plot(par.z/par.units.nm, (par.U_F + par.U_QW + par.dEc * x_ww)/par.units.meV, 'k-', 'LineWidth',2)
        for i = 1 : par.neigs
          plot(par.z/par.units.nm, out.E(i)/par.units.meV + par.N * abs(out.S(:,i)).^2, '-','Color',[1 1 1]*0.5, 'LineWidth',1)  
        end
        plot(par.z/par.units.nm, out.E(out.idx_gnd)/par.units.meV + par.N *  abs(out.S(:,out.idx_gnd)).^2, '-','Color',[1 0 0]*1, 'LineWidth',2)  
        
        xlabel('z (nm)')
        ylabel('U (meV)')
        box on
  
        title(['mean(E_{VS}) = ',num2str(mean_E_VS/par.units.ueV),' ueV | std(E_{VS}) = ',num2str(std_E_VS/par.units.ueV),' ueV'])
  
        %xlim([-2 1]*par.h_QW/par.units.nm)
        xlim([-17, 12])
        ylim([-25, 275])

      %  rice distribution
      plot_counter = plot_counter + 1;      
      subplot(nrows,ncols,1*ncols + iq); hold all;        
        
        x_pdf = linspace(0, mean_E_VS + 5*std_E_VS, 1001)/par.units.ueV;
        pdf   = rician_pdf(x_pdf, nu/par.units.ueV, sigma/par.units.ueV);
        plot(x_pdf, pdf, 'k-','DisplayName','pdf')
        xline(mean_E_VS/par.units.ueV, 'r-','DisplayName','mean')
        box on
        xlabel('energy (ueV)')

        title(['\nu = ',num2str(nu/par.units.ueV),' ueV | \sigma = ',num2str(sigma/par.units.ueV),' ueV | \nu/(2\sigma) = ',num2str(nu/(2*sigma)) ])
        
  
      %  contribution from resonances
      plot_counter = plot_counter + 1;              
      subplot(nrows,ncols,2*ncols + iq); hold all;                
        for in = 1 : length(par.n_range)
          plot(par.n_range(in)*[1 1], [0 1]*abs(Delta_det_n(in))/par.units.ueV, 'b-','LineWidth',5)
        end
        box on
        xlim([-2.5 2.5])
        xlabel('resonance n')
        ylabel('|\Delta_{det,n}|')

      % distribution of Delta in complex plane  
      plot_counter = plot_counter + 1;              
      subplot(nrows,ncols,3*ncols + iq); hold all;                
        plot(real(Delta_det)/par.units.ueV, imag(Delta_det)/par.units.ueV, 'ko')
        phi = linspace(0,2*pi,101);
        %rad = 1/sqrt(2) * sigma/par.units.ueV;
        rad = 1/2 * sigma/par.units.ueV;
        plot( real(Delta_det)/par.units.ueV + rad*cos(phi), imag(Delta_det)/par.units.ueV + rad*sin(phi), 'r-')
        xlim([-1 1]*200)
        ylim([-1 1]*200)
        xlabel('Re(\Delta) (ueV)')
        ylabel('Im(\Delta) (ueV)')   
        xline(0,'k--')
        yline(0,'k--')
        set(gca,'XTick',[-200:100:200])
        set(gca,'YTick',[-200:100:200])
        box on
        axis square


        sigma_gauss = 0.5*sigma/par.units.ueV; % convert sigma from Rice to Normal
        Re_Delta    = linspace(-1,1,501)*200; % in ueV
        Im_Delta    = Re_Delta;
        Re_gaussian = 1/(2*pi*sigma_gauss^2) * exp(-0.5* ((Re_Delta-real(Delta_det/par.units.ueV))/sigma_gauss).^2);
        Im_gaussian = 1/(2*pi*sigma_gauss^2) * exp(-0.5* ((Im_Delta-imag(Delta_det/par.units.ueV))/sigma_gauss).^2);

        surf(Re_Delta,Im_Delta,(Re_gaussian'*Im_gaussian)');
        shading flat
        Ncol = 32;
        cmap = [ones(Ncol,1), linspace(0,1,Ncol)', linspace(0,1,Ncol)'];
        cmap = flipud(cmap);
        %cmap = cmap(Ncol/2+1:end,:);
        colormap(cmap);


      % Fourier trafo
      plot_counter = plot_counter + 1;              
      subplot(nrows,ncols,4*ncols + iq); hold all;                
        U = par.U_QW + par.U_F + x_ww * par.dEc;
        S = U .* out.psi0.^2;
        fftS = fft(S);

        plot( par.k * par.a0/(2*pi), abs(fftS/fftS(1)), 'r-','LineWidth',2)
        set(gca,'XScale','log')
        set(gca,'YScale','log')
        xlabel('k (2\pi/a_0)')
        ylabel('fft(U\psi^2) (a.u.)')   
        box on


        xline(2* par.k1 * par.a0/(2*pi),'m--')
        xline(2* par.k0 * par.a0/(2*pi),'r--')        
        xlim([5E-2,3.5])
        %ylim([3E-6 1E1])
        ylim([5E-3 2E0])
        set(gca,'YTick',[1E-3,1E-2,1E-1,1E0,1E1])


    end
    

  % save
    drawnow
    if save_plot == 1
      height = 800;
      width  = 1600;
      set(fig_obj,'Position',[100 100 width height]);
      exportgraphics(fig_obj,'Fig5cde_map_WW_E_VS_q_vs_x_at_eps_xy_fixed---wavefunctions.pdf','ContentType','vector')
    end

 






end
