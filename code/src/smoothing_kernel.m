function K = smoothing_kernel(k, sigma, model)

  switch(model)
    case 0 % delta
      K = ones(size(k));
    case 1 % Gaussian
      sigma_norm = 2*sqrt(2*log(2));
      sigma = sigma/sigma_norm; 

      K = exp(-0.5*(sigma*k).^2);
    case 2 % Sech
      sigma_norm = 2*acosh(sqrt(2));
      sigma = sigma/sigma_norm; 

      q = (0.5*pi*k*sigma);
      K = q./sinh(q);
      K(q==0) = 1;
    case 3 % Box
      q = sigma*k/2;
      K = sin(q)./q;
      K(q==0) = 1;
    case 4 % Triangle
      q = sigma*k/2;
      K = sin(q)./q;
      K(q==0) = 1;
      K = K.*K;
    case 5 % Lanczos
      q = sigma*k/(2*pi);
      Q = 0.025;
      K = sin(q*pi)./(q*pi) .* sin(q*pi/Q)./(q*pi/Q);
      K(q==0) = 1;      
      K(q>Q)  = 0;
      K(q<-Q) = 0;

    case 6 % Lorentzian
      sigma_norm = 2;
      sigma = sigma/sigma_norm;
      K = exp(-sigma*abs(k));

    otherwise

      error('not implemented')
  end

end

