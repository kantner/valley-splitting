function [c] = cubic_spline_parameters(par)
% cubic spline interpolation for potential in strained Si

% import
  a0       = par.a0;
  %{
  V_0      = par.V_0;
  V_sqrt3  = par.V_sqrt3;
  V_sqrt8  = par.V_sqrt8;
  V_sqrt11 = par.V_sqrt11;
  V_3kF    = par.V_3kF;
  %}

% Kim
  V_0      = -1.113 * par.units.Ry;
  V_sqrt3  = -0.263 * par.units.Ry;
  V_sqrt8  = -0.040 * par.units.Ry;
  V_sqrt11 =  0.033 * par.units.Ry;
  V_3kF    = 0;

  
% extract Fermi level
  EF = -3/2 * V_0;

% Fermi vector
  kF = sqrt(2*par.const.m0*EF)/par.const.hbar;
  
 % cubic spline interpolation on 4 intervals with 4 coefficients per interval
 % s(x) = c1 x^3  + c2 x^2 + c3 x + c4

 % lattice points [in units of a0/(2*pi)]
   x = [0, 2*pi/a0*sqrt(3), 2*pi/a0*sqrt(8), 2*pi/a0*sqrt(11), 3*kF] * a0/(2*pi);

 % matrix and rhs
   A = zeros(16,16);
   b = zeros(16,1);


 %%%%%%%%%%  
 % conditions to match the regular lattice points (5 conditions)
   cond_idx = 1;
   A(cond_idx,1) = x(1)^3;
   A(cond_idx,2) = x(1)^2;
   A(cond_idx,3) = x(1)^1;
   A(cond_idx,4) = 1;
   b(cond_idx)   = V_0;

   cond_idx = 2;
   A(cond_idx,1) = x(2)^3;
   A(cond_idx,2) = x(2)^2;
   A(cond_idx,3) = x(2)^1;
   A(cond_idx,4) = 1;
   b(cond_idx)   = V_sqrt3;
   
   cond_idx = 3;
   A(cond_idx,5) = x(3)^3;
   A(cond_idx,6) = x(3)^2;
   A(cond_idx,7) = x(3)^1;
   A(cond_idx,8) = 1;
   b(cond_idx)   = V_sqrt8;   

   cond_idx = 4;
   A(cond_idx,9)  = x(4)^3;
   A(cond_idx,10) = x(4)^2;
   A(cond_idx,11) = x(4)^1;
   A(cond_idx,12) = 1;
   b(cond_idx)   = V_sqrt11;    

   cond_idx = 5;
   A(cond_idx,13) = x(5)^3;
   A(cond_idx,14) = x(5)^2;
   A(cond_idx,15) = x(5)^1;
   A(cond_idx,16) = 1;
   b(cond_idx)    = V_3kF;    

 % condition to match BC left (V''(0) = 0)
   cond_idx = 6;
   A(cond_idx,1) = 6*x(1);
   A(cond_idx,2) = 2;
   b(cond_idx)   = 0;
   
 % condition to match BC right (V'(3kF) = 0)
   cond_idx = 7;
   A(cond_idx,13) = 3*x(5)^2;
   A(cond_idx,14) = 2*x(5);
   A(cond_idx,15) = 1;
   b(cond_idx)    = 0;   

 % continuity, differentiability and curvature at 1st internal interface x(2)
   cond_idx = 8;
   A(cond_idx,1) = x(2)^3;
   A(cond_idx,2) = x(2)^2;
   A(cond_idx,3) = x(2);
   A(cond_idx,4) = 1;
   A(cond_idx,5) = -x(2)^3;
   A(cond_idx,6) = -x(2)^2;
   A(cond_idx,7) = -x(2);
   A(cond_idx,8) = -1;

   cond_idx = 9;
   A(cond_idx,1) = 3*x(2)^2;
   A(cond_idx,2) = 2*x(2);
   A(cond_idx,3) = 1;
   A(cond_idx,5) = -3*x(2)^2;
   A(cond_idx,6) = -2*x(2);
   A(cond_idx,7) = -1;

   cond_idx = 10;
   A(cond_idx,1) = 2*3*x(2);
   A(cond_idx,2) = 2;
   A(cond_idx,5) = -2*3*x(2);
   A(cond_idx,6) = -2;

 % continuity, differentiability and curvature at 2nd internal interface x(3)
   cond_idx = 11;
   A(cond_idx,5) = x(3)^3;
   A(cond_idx,6) = x(3)^2;
   A(cond_idx,7) = x(3);
   A(cond_idx,8) = 1;
   A(cond_idx,9) = -x(3)^3;
   A(cond_idx,10) = -x(3)^2;
   A(cond_idx,11) = -x(3);
   A(cond_idx,12) = -1;

   cond_idx = 12;
   A(cond_idx,5) = 3*x(3)^2;
   A(cond_idx,6) = 2*x(3);
   A(cond_idx,7) = 1;
   A(cond_idx,9) = -3*x(3)^2;
   A(cond_idx,10) = -2*x(3);
   A(cond_idx,11) = -1;

   cond_idx = 13;
   A(cond_idx,5) = 2*3*x(3);
   A(cond_idx,6) = 2;
   A(cond_idx,9) = -2*3*x(3);
   A(cond_idx,10) = -2;
   
 % continuity, differentiability and curvature at 3rd internal interface x(4)
   cond_idx = 14;
   A(cond_idx,9) = x(4)^3;
   A(cond_idx,10) = x(4)^2;
   A(cond_idx,11) = x(4);
   A(cond_idx,12) = 1;
   A(cond_idx,13) = -x(4)^3;
   A(cond_idx,14) = -x(4)^2;
   A(cond_idx,15) = -x(4);
   A(cond_idx,16) = -1;

   cond_idx = 15;
   A(cond_idx,9) = 3*x(4)^2;
   A(cond_idx,10) = 2*x(4);
   A(cond_idx,11) = 1;
   A(cond_idx,13) = -3*x(4)^2;
   A(cond_idx,14) = -2*x(4);
   A(cond_idx,15) = -1;

   cond_idx = 16;
   A(cond_idx,9) = 2*3*x(4);
   A(cond_idx,10) = 2;
   A(cond_idx,13) = -2*3*x(4);
   A(cond_idx,14) = -2;   


%%%%%%%%%%%
% solve
  c = A\b;

end