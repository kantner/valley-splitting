function C4 = coeff_C4_sym(n_range,tol,par)
% compute bandstructure coefficients C4
% This coefficient is symmetric C_{-n} = C^*_(+n}

  % In the bandstructure calculation, strain is taken into account without
  % explicit smallness assumption. Here it is however more convenient to
  % employ a linear approximation.

  %%%%%%%%%%%%%%%%%%%%%%
  % create list with strain to linear order
  
  % list of integer vectors (no strain)
    g_list = round(par.G_list * par.a0/(2*pi));
  
  % add linear strain
    Id = eye(3);
    for i = 1 : par.N_G
      g_list(:,i) = (Id - par.eps) * g_list(:,i);
    end
  
  %%%%%%%%%%%%%%%%%%%%%%
  % compute bandstructure coefficients
  
  % reduce list of integers to exploit symmetry;
    m_range = unique(abs(n_range));

  % sweep over reduced list

  % allocate memoty
    C4_reduced = zeros(length(m_range),1);  
  
  % scaled shift vector (times 2*pi/a0)
    %g0 = 2*[-par.eps(1,3); - par.eps(2,3); 1-par.eps(3,3)];
    G0 = compute_G0(par);
    g0 = G0 * par.a0/(2*pi);

  % add up coefficients       
    for i1 = 1 : par.N_G
      g1 = g_list(:,i1);

      for i2 = 1 : par.N_G
        g2 = g_list(:,i2);

        for i3 = 1 : par.N_G
          g3 = g_list(:,i3);

          for i4 = 1 : par.N_G
            g4 = g_list(:,i4);
            
            % take sum (using complex conjugation of c_-(G) = c_+^*(-G)
              g = g1 + g2 - g3 - g4;
        
            % run over all integers n
              for im = 1 : length(m_range)
                m = m_range(im);
        
              % check selection rules  
                if abs(g(1) - m*g0(1) ) <= tol
                  if abs(g(2) - m*g0(2) ) <= tol
                    if abs(g(3) - m*g0(3) ) <= tol
                    % add up coefficients  
                      C4_reduced(im) = C4_reduced(im) + par.c(i1)' * par.c(i2)' * par.c(i3) * par.c(i4);
                    end
                  end
                end
              end

          end
        end
      end
    end

  %%%%%%%%%%%%%%%%%%%%%%    
  % build up C_n from C_m
    C4 = zeros(length(n_range),1);  

    for in = 1 : length(n_range)

      n = n_range(in);

      idx = find(abs(n) == m_range);
      if sign(n) < 1
        C4(in) = C4_reduced(idx)'; % complex conjugate
      else
        C4(in) = C4_reduced(idx);
      end

    end



end