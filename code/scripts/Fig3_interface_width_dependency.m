function [] = Fig3_interface_width_dependency(par)
% sweep over interface width


  % exportgraphics
    save_plot = 1; % 0 = off | 1 = on

  % set parameter range
    sigma_range = linspace(0.001,8.5,151) * par.units.ML;

  % DEBUG
    %par.E_cutoff = 4.0 * par.units.Ry;

  % QW
    par.h_QW  = 75 * par.units.ML;
    par.X_barrier = 0.3;
 
  % grid
    par.N     = 2^13;
    par.L     = 5*par.h_QW;
    par.dz    = par.L/par.N;
    par.z     = [0:par.N-1]'*par.dz - par.L/2;

  % k-space grid  
    par.k = [0 : par.N/2-1, -par.N/2:-1]'*2*pi/par.L;

  % set strain
    [eps_QW]     = strain_quantum_well(0.3, par);
    eps_xy_range = [0, 0.1] * 1E-2;

  %%%%%%%%%%%%%%%%%%
  % quantum dot parameters
    par.omega_x = 3.0 * par.units.meV /par.const.hbar;
    par.omega_y = 3.0 * par.units.meV /par.const.hbar;
    par.l_x     = sqrt(par.const.hbar/(par.mass_t*par.omega_x));
    par.l_y     = sqrt(par.const.hbar/(par.mass_t*par.omega_y));

  % set electric field    
    par.F   = 10 * 1E6;

  % potential
    par.U_F = -par.const.e0 * par.F * par.z;    

    
  for ieps = 1 : length(eps_xy_range)
    par.eps = eps_QW;
    par.eps(1,2) = eps_xy_range(ieps);
    par.eps(2,1) = eps_xy_range(ieps);

  % compute conduction band parameters
    [par] = compute_conduction_band_parameters(par.eps, par.xi, par);    



  for is = 1 : length(sigma_range)

    %%%%%%%%%%%%%%%%%  
    % quantum well parameters
      par.sigma_u   = sigma_range(is);
      par.sigma_l   = sigma_range(is);

      
  
    %%%%%%%%%%%%%%%%
    % quantum well indicator
      par.QW_indicator = 0.5 * tanh((par.z+par.h_QW)/par.sigma_l) + 0.5 * tanh(-par.z/par.sigma_u);
  
    % nominal alloy profile    
      par.X_QW         = par.X_barrier * (1 - par.QW_indicator);
      
    % QW potential from nominal (smoothed) profile  
      par.U_QW = par.dEc * par.X_QW;
   
    % wiggle well profile          
      x = zeros(par.N,1);
  
    % compute valley splitting
      compute_derivatives = 0;
      [out] = compute_valley_splitting(x, compute_derivatives, par);
      
      
      mean_E_VS = out.M;
      var_E_VS  = out.V;
      nu        = out.nu;
      sigma     = out.sigma;
      Delta_det = out.Delta_det;

      Delta_det_n = out.Delta_det_n;
          
  
      map.mean_E_VS(ieps,is) = mean_E_VS;
      map.std_E_VS(ieps,is)  = sqrt(var_E_VS);
      map.nu(ieps,is)        = nu;
      map.sigma(ieps,is)     = sigma; 
      map.Delta_det(ieps,is) = Delta_det; 
      map.Delta_det_n(ieps,is,:) = Delta_det_n; 
     
  end
  end

 %% compute confidence intervals
    upper_bound = 0.75;
    lower_bound = 0.25;   

    for ieps = 1 : 2
      for is = 1 : length(sigma_range)
        
        nu    = map.nu(ieps,is);
        sigma = map.sigma(ieps,is);
      
        map.sigma_lower(ieps,is) = par.units.ueV * compute_rice_percentile(lower_bound, nu/par.units.ueV, sigma/par.units.ueV);
        map.sigma_upper(ieps,is) = par.units.ueV * compute_rice_percentile(upper_bound, nu/par.units.ueV, sigma/par.units.ueV);


      end
    end

 %% %%%%%%%%%%%
  % plot
    fig_obj = figure(60001); clf; hold all;
    sgtitle('interface width dependency')

    ymin = -20;
    ymax = 320;

    length_scale = par.units.ML;
    energy_scale = par.units.ueV;

    for ieps = 1 : 2
      subplot(2,2,ieps); hold all;
      plot(sigma_range/length_scale,map.mean_E_VS(ieps,:)/energy_scale,'k-','LineWidth',2,'DisplayName','mean(E_{VS})')
      plot(sigma_range/length_scale,map.std_E_VS(ieps,:)/energy_scale,'r-','LineWidth',2,'DisplayName','std(E_{VS})')

      %plot(sigma_range/length_scale,map.mean_E_VS(ieps,:)/energy_scale,'k-','LineWidth',2,'DisplayName','mean(E_{VS})')

      plot(sigma_range/length_scale,map.nu(ieps,:)/energy_scale,'c-','LineWidth',2,'DisplayName','\nu')
      plot(sigma_range/length_scale,map.sigma(ieps,:)/energy_scale,'m-','LineWidth',2,'DisplayName','\sigma')
      

      patch_obj = patch([sigma_range, flip(sigma_range)]/length_scale, [map.sigma_lower(ieps,:), flip(map.sigma_upper(ieps,:))]/energy_scale,'r');
      set(patch_obj,'FaceColor',[1 0 0]*1,'FaceAlpha',0.25,'EdgeColor','none','DisplayName','(25%-75%) quantile')

      %box on
      xlabel('interface width (ML)')
      ylabel('energy (ueV)')
      
      legend('Location','best')
      ylim([ymin ymax])
      xlim([min(sigma_range) max(sigma_range)]/length_scale)
      %set(gca,'XTick',[0:0.25:1.6])
      set(gca,'XTick',[0:1:16])

      
  
    %%%%%%%%%%%%  
    % find cross-over point for deterministic enhancement
      sigma_crit = interp1(map.nu(ieps,:)./(0.3507 * 2*map.sigma(ieps,:)),sigma_range,1,'pchip');  
      xline(sigma_crit/length_scale,'k--','DisplayName','det. enhancement')
    

      

    %%%%%%  
    % second axis (ML = monolayer) 
    %{
      ML = par.a0/4;
      hax(1) = gca;
      hax(2) = axes('Position', hax(1).Position, 'XAxisLocation', 'top', 'YAxisLocation','right','color','none');
      hold(hax(2),'on')
      
      ylim([ymin ymax])
      xlim([min(sigma_range) max(sigma_range)]/ML)
      
    % plot single line 
      plot(hax(2), sigma_range/ML, map.mean_E_VS(ieps,:)/par.units.ueV)
      

      set(gca,'XTick',[0:16])
      xlabel('interface width (ML)')
    %}

      title(['\epsilon_{x,y} = ',num2str(eps_xy_range(ieps)*100),'%'])
      box on
      ylim([ymin ymax])


      subplot(2,2,ieps+2); hold all;
      cmap = lines(length(par.n_range));
      
      plot(sigma_range/length_scale,abs(map.Delta_det(ieps,:))/energy_scale,'-.','Color',[0 0 0],'LineWidth',4,'DisplayName',['|\Delta_{det}|'])
      for in = 1 : length(par.n_range)
        if sum(abs(map.Delta_det_n(ieps,:,in))) > 0
        plot(sigma_range/length_scale,abs(map.Delta_det_n(ieps,:,in))/energy_scale,'-','Color',cmap(in,:),'LineWidth',2,'DisplayName',['|\Delta_{det,n}| (n=',num2str(par.n_range(in)),')'])
        end        
      end
      
      legend()
      box on;
      set(gca,'YScale','log')
      ylim([0.3E-1, 8E2])      
      xlim([min(sigma_range) max(sigma_range)]/length_scale)
      xlabel('interface width (ML)')
      ylabel('energy (ueV)')
      set(gca,'XTick',[0:1:16])

  
    end
    %%%%%%%% 
    % save
      drawnow
      if save_plot == 1
        height = 800;
        width  = 1600;
        set(fig_obj,'Position',[100 100 width height]);
        exportgraphics(fig_obj,'Fig3_interface_width_dependency.pdf','ContentType','vector')
      end  
  






end
