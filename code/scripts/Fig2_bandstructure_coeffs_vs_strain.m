function [] = Fig2_bandstructure_coeffs_vs_strain(par)

  %%%%%%%%%%%%%
  % energy cutoff
    par.E_cutoff = 8 * par.units.Ry; % reduce for fast debugging


  % exportgraphics
    save_plot = 1; % 0 = off | 1 = on


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % shear strain range
    eps_range = linspace(0,0.5,21) * 1E-2;
    
  % allocate memory
    k0_list = zeros(1, length(eps_range));
    C2_list = zeros(length(par.n_range), length(eps_range));
    C4_list = zeros(length(par.n_range), length(eps_range));

  % quantum well strain
    [eps_QW] = strain_quantum_well(0.3, par);
     
  
  % compute   
    for ieps = 1 : length(eps_range)

      % set up strain tensor
        eps_xy = eps_range(ieps);

        par.eps      = eps_QW;
        par.eps(1,2) = eps_xy;
        par.eps(2,1) = eps_xy;

      % compute conduction band parameters
        [par] = compute_conduction_band_parameters(par.eps, par.xi, par);  

      % store k0
        k0_list(ieps)   = par.k0_vec(3); 
        C2_list(:,ieps) = par.C2;
        C4_list(:,ieps) = par.C4;
    end    

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % bandstructure coefficients C2
    fig_obj = figure(90004); clf; hold all;

      sgtitle('Figure 2(a)')

      for i = 1 : 2
      subplot(1,2,i); hold all;

      cmap = lines(length(par.n_range));
      for in = 1 : length(par.n_range)
        if sum(abs(C2_list(in,:))) > 0
        plot(eps_range*100, C2_list(in,:),'-','Color',cmap(in,:),'LineWidth',2,'DisplayName',['C^{(2)}_{',num2str(par.n_range(in)),'}'])
        end
      end

      plot(eps_range*100, -58.9*eps_range,'m--','Color',cmap(in,:),'LineWidth',2,'DisplayName',['fit'])


      legend('Location','eastoutside')
      box on
      xlabel('\epsilon_{x,y} (%)')
      ylabel('band structure coefficients C^{(2)}')
        if i == 1
          set(gca,'XScale','linear')
          set(gca,'YScale','linear')
          ylim([-0.005 0.001])
        else
          set(gca,'XScale','linear')
          set(gca,'YScale','linear')
          ylim([-0.31 -0.04])
        end
      end


    % save
      drawnow
      if save_plot == 1
        height = 800;
        width  = 1600;
        set(fig_obj,'Position',[100 100 width height]);
        exportgraphics(fig_obj,'Fig2a_C2_vs_eps_xy.pdf','ContentType','vector')
      end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % bandstructure coefficients C4m
    fig_obj = figure(90005); clf; hold all;

      sgtitle('Figure 2(b)')

      for i = 1 : 2
      subplot(1,2,i); hold all;

      cmap = lines(length(par.n_range));
      for in = 1 : length(par.n_range)
        if sum(abs(C4_list(in,:))) > 0  && par.n_range(in) >= 0 % only positive because of symmetry     
        plot(eps_range*100, C4_list(in,:),'-','Color',cmap(in,:),'LineWidth',2,'DisplayName',['C^{(4,-)}_{\pm',num2str(par.n_range(in)),'}'])
        end
      end
      legend('Location','eastoutside')
      box on
      xlabel('\epsilon_{x,y} (%)')
      ylabel('band structure coefficients C^{(4,-)}')
        if i == 1
          set(gca,'XScale','linear')
          set(gca,'YScale','linear')
          ylim([-0.01 0.145])  
        elseif i == 2
          set(gca,'XScale','linear')
          set(gca,'YScale','linear')
          ylim([1.3 1.4])  
        end
      end
      
    % save
      drawnow
      if save_plot == 1
        height = 800;
        width  = 1600;
        set(fig_obj,'Position',[100 100 width height]);        
        exportgraphics(fig_obj,'Fig2b_C4m_vs_eps_xy.pdf','ContentType','vector')
      end    

end