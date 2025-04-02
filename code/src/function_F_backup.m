function [F] = function_F(K1,K2,R)
% function entering the non-local pseudopotential

  X1 = K1 * R;
  X2 = K2 * R;


if K1 == 0 && K2 == 0
  F = 1/3 * R^3 ;
elseif K1 == 0 && K2 > 0
  X = K2 * R;
  F = R^3 * spherical_besselj(1,X)/X;
elseif K2 == 0 && K1 > 0
  X = K1 * R;
  F = R^3 * spherical_besselj(1,X)/X;

elseif K1 == K2 && K1 > 0
%    tmp =  spherical_besselj(0, K1*R);
%    F= 0.5*R^3 *(tmp^2 -  spherical_besselj(-1, K1*R) *  spherical_besselj(+1, K1*R));



    X = K1*R;
    F= 0.5*R^3 *(1 - sin(2*X)/(2*X))/(X*X);
  else
    F= R^2/(K1^2 - K2^2) *(K1 * spherical_besselj(+1, K1*R) * spherical_besselj(0, K2*R) - K2 * spherical_besselj(+1, K2*R)* spherical_besselj(0, K1*R));
end



end