function [F] = function_F(K1,K2,R)
% function entering the non-local pseudopotential

  X1 = K1 * R;
  X2 = K2 * R;


 % cases
 % 1: X1 \approx X2 
 % 2: X1 \neq X2 

 eps_threshold = 1E-3;

 X_small_threshold = 1E-5;

 if abs(X1-X2) < eps_threshold
 
   X   = (X1+X2)/2; % mean
   eps = X1-X2;   % difference


   if X < X_small_threshold
   % both X1 = X2 are very small
    f = 1/3 - X*X/15 + eps*eps *(-1/60 + X*X/(630));
   else
   % both X1 = X2 are equal but not small
    Xpow3 = X*X*X;
    Xpow5 = Xpow3 * X*X;
    sin2X = sin(2*X);
    f  = (X - 0.5*sin2X)/(2*Xpow3) + eps*eps * (3*X-2*Xpow3-3*0.5*sin2X)/(24*Xpow5); %2nd order Taylor expansion
   end
 
 else

   denom = (X1-X2)*(X1+X2);
   f = (X1 * spherical_besselj(1,X1) * spherical_besselj(0,X2) - X2*spherical_besselj(1,X2)*spherical_besselj(0,X1))/denom;

 end

 F = R^3 * f;

end
