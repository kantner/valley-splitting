function [V] = V_atomic(G, par)


  switch(par.V_atomic_model)
    case 1
      V = V_atomic_Fischetti_Laux(G, par);
    case 2
      V = V_atomic_Kim_Fischetti(G, par);
    case 3
      V = V_atomic_Fischetti_Higman(G, par);
    case 4
      V = V_atomic_Friedel(G, par);
    otherwise
      error('invalid option')
  end

end