function [] = Fig4_wiggle_well_line_plots(par)
% make plot like Joynt/Feng (2022), Fig 3

  % exportgraphics
    save_plot = 1; % 0 = off | 1 = on

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % strain tensor for QW
    [eps_QW] = strain_quantum_well(0.3, par);

  % do not compute derivates of valley-splitting functions
    compute_derivatives = 0;

  %%%%%%%%%%%%%%%%%  
  % quantum well parameters
    par.h_QW      = 75  * par.units.ML;
    par.X_barrier = 0.3;
    par.sigma_u   = 0.5 * par.units.nm;
    par.sigma_l   = 0.5 * par.units.nm;

  % electric field (in V/m)
    par.F       = 5.0 * 1E6; % considered in Lima 2023

  %%%%%%%%%%%%%%%%
  % grid
    par.N  = 2^12;
    par.L  = 5*par.h_QW;
    par.dz = par.L/par.N;
    par.z  = [0:par.N-1]'*par.dz - par.L/2;

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
    

  % k-space grid  
    par.k = [0 : par.N/2-1, -par.N/2:-1]'*2*pi/par.L;


  % create fig object
    fig_obj = figure(700001); clf; hold all;    
    nrows = 3;
    ncols = 1;


  %%%%%%%%%%%
  % FIRST PLOT
  % no shear strain vs. q, different x
    
    eps_xy   = 0 * par.units.percent;
    q_range  = linspace(0,2,301) * 2*pi/par.a0;
    x_range  = [0.01 0.05 0.1 0.15];

  % set QW strain + fixed shear strain 
    eps      = eps_QW;
    eps(1,2) = eps_xy;
    eps(2,1) = eps_xy;
    
  % compute conduction band parameters    
    [par] = compute_conduction_band_parameters(eps, par);

  % sweep  
    for ix = 1 : length(x_range)
      for iq = 1 : length(q_range)

      % wiggle well profile          
        x_ww = x_wiggle_well(x_range(ix), q_range(iq), 0, par);

      % compute valley splitting
        [out] = compute_valley_splitting(x_ww, eps, compute_derivatives, par);
    
        mean_E_VS = out.M;
        var_E_VS  = out.V;
        nu        = out.nu;
        sigma     = out.sigma;
        

        map(1).mean_E_VS(iq, ix) = mean_E_VS;
        map(1).std_E_VS(iq, ix)  = sqrt(var_E_VS);
        map(1).nu(iq, ix)        = nu;
        map(1).sigma(iq, ix)     = sigma;        

      end

    end



  %%%%%%%%%%%
  % SECOND PLOT
  % small shear strain vs. q, different x
    eps_xy   = 0.1 * par.units.percent;
    
  % set QW strain + fixed shear strain 
    eps      = eps_QW;
    eps(1,2) = eps_xy;
    eps(2,1) = eps_xy;
    
  % compute conduction band parameters
    [par] = compute_conduction_band_parameters(eps, par);

  % sweep  
    for ix = 1 : length(x_range)

      for iq = 1 : length(q_range)

      % wiggle well profile          
        x_ww = x_wiggle_well(x_range(ix), q_range(iq), 0, par);

      % compute valley splitting
        [out] = compute_valley_splitting(x_ww, eps, compute_derivatives, par);
    
        mean_E_VS = out.M;
        var_E_VS  = out.V;
        nu        = out.nu;
        sigma     = out.sigma;
        

        map(2).mean_E_VS(iq, ix) = mean_E_VS;
        map(2).std_E_VS(iq, ix)  = sqrt(var_E_VS);
        map(2).nu(iq, ix)        = nu;
        map(2).sigma(iq, ix)     = sigma;        

      end

    end


            
  %%%%%%%%%%%
  % THIRD PLOT
  % fixed amplitude, different strains vs. q
    
    %eps_xy_range = [0 0.05 0.1 0.2] * par.units.percent;
    eps_xy_range = [0 0.1 0.2 0.5] * par.units.percent;
    x_max        = 0.05;

    for ieps = 1 : length(eps_xy_range)

      % set QW strain + fixed shear strain 
        eps_xy = eps_xy_range(ieps);
    
        eps      = eps_QW;
        eps(1,2) = eps_xy;
        eps(2,1) = eps_xy;

      % compute conduction band parameters
        [par] = compute_conduction_band_parameters(eps, par);

      for iq = 1 : length(q_range)

      % wiggle well profile          
        x_ww = x_wiggle_well(x_max, q_range(iq), 0, par);

      % compute valley splitting
        [out] = compute_valley_splitting(x_ww, eps, compute_derivatives, par);
    
        mean_E_VS = out.M;
        var_E_VS  = out.V;
        nu        = out.nu;
        sigma     = out.sigma;
        

        map(3).mean_E_VS(iq, ieps) = mean_E_VS;
        map(3).std_E_VS(iq, ieps)  = sqrt(var_E_VS);
        map(3).nu(iq, ieps)        = nu;
        map(3).sigma(iq, ieps)     = sigma;        

      end

    end

  %% make plots
  %%%%%%%%%%%%%%%%%%%%% 
    fig_obj = figure(700001); clf; hold all;    


  % plot 1    
    subplot(nrows,ncols,1); hold all;
    cmap = lines(length(x_range));
    for ix = 1 : length(x_range)
      plot(q_range*par.a0/(2*pi), map(1).mean_E_VS(:,ix)/par.units.meV, '-','Color',cmap(ix,:),'LineWidth',2,'DisplayName',['mean(E_{VS}) at x = ',num2str(x_range(ix))])
      %plot(q_range*par.a0/(2*pi), map(1).std_E_VS(:,ix)/par.units.meV, '--','Color',cmap(ix,:),'LineWidth',2,'DisplayName',['std(E_{VS})  at x = ',num2str(x_range(ix))])
      %plot(q_range*par.a0/(2*pi), map(1).nu(:,ix)/par.units.meV, '--','Color',cmap(ix,:),'LineWidth',2,'DisplayName',['std(E_{VS})  at x = ',num2str(x_range(ix))])
      %plot(q_range*par.a0/(2*pi), map(1).sigma(:,ix)/par.units.meV, '--','Color',cmap(ix,:),'LineWidth',2,'DisplayName',['std(E_{VS})  at x = ',num2str(x_range(ix))])
    end
    xlabel('q (2\pi/a_0)')    
    ylabel('energy (meV)')      
    title(['no shear strain, various Ge amplitudes'])
    box on
    legend('Location','eastoutside')
    ylim([-1.0 9.0])
    xline(2*par.k0 * par.a0/(2*pi),'m--','Linewidth',2,'DisplayName','2k_0')
    xline(2*par.k1 * par.a0/(2*pi),'c--','Linewidth',2,'DisplayName','2k_1')



  % plot 2
    fig_obj = figure(700001); hold all;    
    subplot(nrows,ncols,2); hold all;
    cmap = lines(length(x_range));
    for ix = 1 : length(x_range)
      plot(q_range*par.a0/(2*pi), map(2).mean_E_VS(:,ix)/par.units.meV, '-','Color',cmap(ix,:),'LineWidth',2,'DisplayName',['mean(E_{VS}) at x = ',num2str(x_range(ix))])
      %plot(q_range*par.a0/(2*pi), map(2).std_E_VS(:,ix)/par.units.meV, '--','Color',cmap(ix,:),'LineWidth',2,'DisplayName',['std(E_{VS})  at x = ',num2str(x_range(ix))])
    end
    xlabel('q (2\pi/a_0)')    
    ylabel('energy (meV)')      
    title(['fixed shear strain 0.1% , various Ge amplitudes'])
    box on
    legend('Location','eastoutside')    
    ylim([-1.0 9.0])

    xline(2*par.k0 * par.a0/(2*pi),'m--','Linewidth',2,'DisplayName','2k_0')
    xline(2*par.k1 * par.a0/(2*pi),'c--','Linewidth',2,'DisplayName','2k_1')

  % plot 3
    fig_obj = figure(700001); hold all;    
    subplot(nrows,ncols,3); hold all;
    cmap = lines(length(eps_xy_range));
    for ieps = 1 : length(eps_xy_range)
      plot(q_range*par.a0/(2*pi), map(3).mean_E_VS(:,ieps)/par.units.meV, '-','Color',cmap(ieps,:),'LineWidth',2,'DisplayName',['mean(E_{VS}) at \epsilon_{x,y} = ',num2str(eps_xy_range(ieps)*100),'%'])
      plot(q_range*par.a0/(2*pi), map(3).std_E_VS(:,ieps)/par.units.meV, '--','Color',cmap(ieps,:),'LineWidth',2,'DisplayName',['std(E_{VS})  at \epsilon_{x,y} = ',num2str(eps_xy_range(ieps)*100),'%'])

     % plot(q_range*par.a0/(2*pi), map.mean_E_VS(:,ieps)/par.units.meV, '-','Color',cmap(ieps,:),'LineWidth',2,'DisplayName',['std(E_{VS})  at \epsilon_{x,y} = ',num2str(eps_xy_range(ieps)*100),'%'])
     % plot(q_range*par.a0/(2*pi), (map.mean_E_VS(:,ieps) + map.std_E_VS(:,ieps))/par.units.meV, '--','Color',cmap(ieps,:),'LineWidth',2,'DisplayName',['std(E_{VS})  at \epsilon_{x,y} = ',num2str(eps_xy_range(ieps)*100),'%'])
     % plot(q_range*par.a0/(2*pi), (map.mean_E_VS(:,ieps) - map.std_E_VS(:,ieps))/par.units.meV, '--','Color',cmap(ieps,:),'LineWidth',2,'DisplayName',['std(E_{VS})  at \epsilon_{x,y} = ',num2str(eps_xy_range(ieps)*100),'%'])


     %{
      patch_obj = patch([q_range,fliplr(q_range)]*par.a0/(2*pi),...
                        ([map.mean_E_VS(:,ieps) + map.std_E_VS(:,ieps);flipud(map.mean_E_VS(:,ieps) - map.std_E_VS(:,ieps))])/par.units.meV, 'r');

      set(patch_obj, 'EdgeColor','none','FaceColor',cmap(ieps,:),'FaceAlpha',0.2)
     %}
    end
    xline(2*par.k0 * par.a0/(2*pi),'m--','Linewidth',2,'DisplayName','2k_0')
    xline(2*par.k1 * par.a0/(2*pi),'c--','Linewidth',2,'DisplayName','2k_1')
    
    xlabel('q (2\pi/a_0)')    
    ylabel('energy (meV)')      
    title(['various shear strains, fixed Ge amplitude x = ',num2str(x_max)])
    box on
    legend('Location','eastoutside')
    ylim([-0.25 3.25])
    

  % save
    drawnow
    if save_plot == 1
      height = 800;
      width  = 1600;
      set(fig_obj,'Position',[100 100 width height]);
      exportgraphics(fig_obj,'Fig4_WW_line_scan_q_fixed_x_different_eps_xy.pdf','ContentType','vector')
    end




end




