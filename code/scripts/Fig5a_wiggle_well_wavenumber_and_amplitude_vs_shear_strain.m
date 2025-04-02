function [] = Fig5a_wiggle_well_wavenumber_and_amplitude_vs_shear_strain(par)

  %%%%%%%%%%%%%
  % energy cutoff
    par.E_cutoff = 8 * par.units.Ry; % reduce for fast debugging


  % exportgraphics
    save_plot = 1; % 0 = off | 1 = on


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % shear strain range
    eps_range = linspace(0,0.2,101) * 1E-2;
    
  % wave number and amplitude range   
    x_range   = linspace(0,0.2,251);           % Fig 2
    
 
  % fixed k1
    k1 = 0.16* 2*pi/par.a0;
    
  % allocate memory
    k0_list  = zeros(1, length(eps_range));
    C2_list  = zeros(length(par.n_range), length(eps_range));
    C4m_list = zeros(length(par.n_range), length(eps_range));

  % quantum well strain
    [eps_QW] = strain_quantum_well(0.3, par);

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
    par.N  = 2^11;
    par.L  = 5*par.h_QW;
    par.dz = par.L/par.N;
    par.z  = [0:par.N-1]'*par.dz - par.L/2;

  % k-space grid  
    par.k = [0 : par.N/2-1, -par.N/2:-1]'*2*pi/par.L;

  %%%%%%%%%%%%%%%%
  % quantum well indicator
    par.QW_indicator = 0.5 * tanh((par.z+par.h_QW)/par.sigma_l) + 0.5 * tanh(-par.z/par.sigma_u);

  % nominal alloy profile    
    par.X_QW         = par.X_barrier * (1 - par.QW_indicator);
    
  % QW potential from nominal (smoothed) profile  
    par.U_QW = par.dEc * par.X_QW;

  % electric field potential
    par.U_F  = -par.const.e0 * par.F * par.z;        
  
  % compute map  
    for ieps = 1 : length(eps_range)

      % set up strain tensor
        eps_xy = eps_range(ieps);

        par.eps      = eps_QW;
        par.eps(1,2) = eps_xy;
        par.eps(2,1) = eps_xy;

      % compute conduction band parameters
        [par] = compute_conduction_band_parameters(par.eps, par.xi, par);  

      % store k0
        k0_list(ieps)    = par.k0_vec(3); 
        C2_list(:,ieps)  = par.C2;
        C4m_list(:,ieps) = par.C4;

      for ix = 1 : length(x_range)
      
      % wiggle well profile          
        x_ww = x_wiggle_well(x_range(ix), 2*k1, 0, par);

      % compute valley splitting
        compute_derivatives = 0;
        [out] = compute_valley_splitting(x_ww, compute_derivatives, par);
    
        mean_E_VS = out.M;        
        std_E_VS  = sqrt(out.V);
        nu        = out.nu;
        sigma     = out.sigma;
        
        map.mean_E_VS(ix,ieps) = mean_E_VS;
        map.std_E_VS(ix,ieps)  = std_E_VS;
        map.nu(ix,ieps)        = nu;
        map.sigma(ix,ieps)     = sigma;

      end

    end


%%  %%%%%%%%%%%%%%%%%%%%%%%
    % plot    
      x_axis_range = x_range;
      x_axis_label = 'x (Ge amplitude)';
      plot_title   = ['fixed q = 2\timesk_1 = ',num2str(2*k1 * par.a0/(2*pi)),' \times 2\pi/a_0'];
      file_name    = 'Fig5a_map_WW_E_VS_x_vs_eps_at_q_fixed.pdf';
      x_axis_min   = 0;
      x_axis_max   = 0.2;
      c_axis_min   = 0;
      c_axis_max   = 3.5;
      %cmap         = divergingColormap1(25);
      cmap         = flipud(lajolla(14));

      y_axis_range = eps_range * 100;
      y_axis_label = 'shear strain \epsilon_{x,y} (%)';
      y_axis_min   = 0;
      y_axis_max   = 0.15;
        
      fig_obj = figure(80000); clf; hold all;
        nrows = 2;
        ncols = 3;
        sgtitle(plot_title)


  
      subplot(nrows,ncols,1); hold all;
        surf(x_axis_range, y_axis_range, (map.mean_E_VS/par.units.meV)')
        xlabel(x_axis_label)
        ylabel(y_axis_label)
        zlabel('mean E_{VS} (meV)')      
        title(['mean(E_{VS})'])
        colorbar
        box on
        shading interp
        colormap(cmap)
        xlim([x_axis_min x_axis_max])
        ylim([y_axis_min y_axis_max])
        clim([c_axis_min c_axis_max])
        colormap(cmap)  
  
      subplot(nrows,ncols,2); hold all;
        surf(x_axis_range, y_axis_range, (map.std_E_VS/par.units.meV)')
        xlabel(x_axis_label)
        ylabel(y_axis_label)
        zlabel('std E_{VS} (meV)')
        title(['std(E_{VS})'])
        colorbar
        box on
        shading interp
        xlim([x_axis_min x_axis_max])
        ylim([y_axis_min y_axis_max])
        clim([c_axis_min c_axis_max])
        colormap(cmap)  
  
      subplot(nrows,ncols,3); hold all;
        surf(x_axis_range, y_axis_range, (map.std_E_VS./map.mean_E_VS)')
        xlabel(x_axis_label)
        ylabel(y_axis_label)
        zlabel('mean E_{VS} (meV)')      
        title(['ratio std(E_{VS})/mean(E_{VS})'])
        colorbar
        box on
        shading interp
        xlim([x_axis_min x_axis_max])
        ylim([y_axis_min y_axis_max])
        clim([c_axis_min c_axis_max])
        colormap(cmap)  
    
      subplot(nrows,ncols,4); hold all;
        surf(x_axis_range, y_axis_range, (map.nu/par.units.meV)')
        xlabel(x_axis_label)
        ylabel(y_axis_label)        
        zlabel('mean E_{VS} (meV)')
        title(['Rice parameter \nu'])
        colorbar    
        box on
        shading interp
        xlim([x_axis_min x_axis_max])
        ylim([y_axis_min y_axis_max])
        clim([c_axis_min c_axis_max])
        colormap(cmap)  
  
      subplot(nrows,ncols,5); hold all;
        surf(x_axis_range, y_axis_range, (map.sigma/par.units.meV)')
        xlabel(x_axis_label)
        ylabel(y_axis_label)
        zlabel('mean E_{VS} (meV)')
        title(['Rice parameter \sigma'])
        colorbar   
        box on
        shading interp
        colormap(cmap)
        xlim([x_axis_min x_axis_max])
        ylim([y_axis_min y_axis_max])
        clim([c_axis_min c_axis_max])
        colormap(cmap)  
  
      subplot(nrows,ncols,6); hold all;
        surf(x_axis_range, y_axis_range, (map.nu./(2*map.sigma))')
        xlabel(x_axis_label)
        ylabel(y_axis_label)
        zlabel('\nu/(2\sigma)')
        title(['ratio of Rice parameters \nu/(2\sigma)'])
        colorbar   
        box on      
        shading interp
        xlim([x_axis_min x_axis_max])
        ylim([y_axis_min y_axis_max])
        clim([c_axis_min c_axis_max])
        colormap(cmap)  
  
  
      %%%%%%%%%%%%%%%%%%%%%%%%
      % contour
        contour_levels = [0.0 1.0]; % we need to add a second contour level, because somehow it does not do it for a single contour level alone (add 0 as dummy)
        M = contour(x_axis_range, y_axis_range, (map.nu./(0.3507 * 2*map.sigma))',contour_levels);
         
        
      % total array size
        n = size(M,2);
    
        counter = 1;
    
      % plot all contour lines
        while counter <= n
    
          val   = M(1,counter);
          n_seg = M(2,counter);
          x = M(1,counter+1: counter+n_seg);
          y = M(2,counter+1: counter+n_seg);
      
          if val == 1 % plot only the val=1 contour line
            for i = 1 : nrows*ncols
              subplot(nrows,ncols,i); hold all;  
              plot3(x,y,val * ones(n_seg,1),'k-','LineWidth',2)
            end
          end
    
          counter = counter+n_seg+1;
        end


    
    %%%%%%%%%%%%%%
    % add lines
    %{
        for i = 1 : nrows*ncols
          subplot(nrows,ncols,i); hold all;  
          z_dummy = 100 * ones(size(y_axis_range));
          plot3(2*k0_list * par.a0/(2*pi)      , y_axis_range, z_dummy,'w--','LineWidth',2)
          plot3(2*(1 - k0_list * par.a0/(2*pi)), y_axis_range, z_dummy,'w--','LineWidth',2)  
          plot3(1*(1 - k0_list * par.a0/(2*pi)), y_axis_range, z_dummy,'m--','LineWidth',2)  
        end
    %}
      

    %%%%%%%% 
    % save
      drawnow
      if save_plot == 1
        height = 800;
        width  = 1600;
        set(fig_obj,'Position',[100 100 width height]);
        exportgraphics(fig_obj,file_name,'ContentType','vector')
      end





end