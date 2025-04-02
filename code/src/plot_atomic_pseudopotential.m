function [] = plot_atomic_pseudopotential(par)
  
    G_list = sort([linspace(0,25,1001),sqrt(3),sqrt(8),sqrt(11),sqrt(12),sqrt(16)])* 2*pi/par.a0;    
    V_a   = V_atomic_Fischetti_Laux(G_list, par);
    V_a2  = V_atomic_Kim_Fischetti(G_list, par);
    V_a3  = V_atomic_Fischetti_Higman(G_list, par);    
    V_a4  = V_atomic_Friedel(G_list, par);


    figure(24354);clf;hold all;
      plot(G_list*par.a0/(2*pi), V_a/par.units.Ry,'ko-','DisplayName','Fischetti/Laux')
      plot(G_list*par.a0/(2*pi), V_a2/par.units.Ry,'bx-','DisplayName','Kim/Fischetti')
      plot(G_list*par.a0/(2*pi), V_a3/par.units.Ry,'ms-','DisplayName','Fischetti/Higman')
      plot(G_list*par.a0/(2*pi), V_a4/par.units.Ry,'gd-','DisplayName','Friedel')
      xline(sqrt(3),'r-')
      yline(-0.263,'r-')

      xline(sqrt(8),'g-')
      yline(-0.040,'g-')

      xline(sqrt(11),'b-')
      yline(+0.033,'b-')

      xline(sqrt(12),'c-')
      xline(sqrt(16),'c-')
      box on
      xlabel('lattice vector modulus G (2\pi/a_0)')
      ylabel('atomic potential V_a (Ry)')
      title('atomic potential (Fischetti/Laux)')
      ylim([-1 0.2])
      xlim([0 5])
      legend()
      drawnow

    % atomic form factors  
    %{
      V_atomic(sqrt(3)*2*pi/par.a0, par)/par.units.Ry
      V_atomic(sqrt(8)*2*pi/par.a0, par)/par.units.Ry
      V_atomic(sqrt(11)*2*pi/par.a0, par)/par.units.Ry
  %}

end
