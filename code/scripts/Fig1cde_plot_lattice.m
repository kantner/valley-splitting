   
 % set strain
   eps = zeros(3,3);
   xi = 0.53;

 % plot lattice for different strains  
   figure(40001);clf;hold all;
      subplot(1,3,1); hold all;
        eps = zeros(3,3);
        [R_list] = plot_diamond_lattice(eps,xi)

      subplot(1,3,2); hold all;
        eps = zeros(3,3);
        eps(1,1) = 0.15;
        eps(2,2) = eps(1,1);
        eps(3,3) = -0.15;
        [R_list] = plot_diamond_lattice(eps,xi)

      subplot(1,3,3); hold all;
        eps(1,2) = -0.2;
        eps(2,1) = eps(1,2);
        [R_list] = plot_diamond_lattice(eps,xi)


      drawnow
  %exportgraphics(gcf,'strained_lattice.pdf','ContentType','vector')
  exportgraphics(gcf,'Fig1b_strained_lattice.png','ContentType','image','Resolution',1200,'BackgroundColor','white')

 %% change to top view and export again
    for i = 1 : 3
      subplot(1,3,i)
      view([0 90])
    end

    drawnow

    exportgraphics(gcf,'Fig1b_strained_lattice_top_view.png','ContentType','image','Resolution',1200,'BackgroundColor','white')


 %% plot primitive unit cells
    figure(40002);clf; hold all;

      subplot(1,3,1); hold all;
        eps = zeros(3,3);
        plot_primitive_unit_cell(eps,xi)

      subplot(1,3,2); hold all;
        eps = zeros(3,3);
        eps(1,1) = 0.3;
        eps(2,2) = eps(1,1);
        eps(3,3) = -0.5;
        plot_primitive_unit_cell(eps,xi)

      subplot(1,3,3); hold all;
        eps(1,2) = -0.2;
        eps(2,1) = eps(1,2);
        plot_primitive_unit_cell(eps,xi)

