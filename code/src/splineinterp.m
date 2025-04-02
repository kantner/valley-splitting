
% cubic spline interpolation

physical_constants;
Ry = rydberg;
m0 = electronMass;
hbar = reducedPlanckConstant;

V_0      = -1.113 * Ry;
V_sqrt3  = -0.263 * Ry;
V_sqrt8  = -0.040 * Ry;
V_sqrt11 = +0.033 * Ry;
V_3kF    = 0;

nm = 1E-9;
a0 = 0.541*nm;

EF = -3/2 * V_0;
kF = sqrt(2*m0*EF)/hbar;
  
 % cubic spline interpolation on 4 intervals with 4 coefficients per interval
 % s(x) = c1 x^3  + c2 x^2 + c3 x + c4

 % lattice points
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


%%%%%%%%%%%
% plot
  figure(1); clf; hold all;
  for i = 1 : 4
  q_range = linspace(x(i), x(i+1), 101); %* 2*pi/a0;
  Vq = c((i-1)*4+1)  * q_range.^3 + c((i-1)*4+2) * q_range.^2 + c((i-1)*4+3) * q_range + c((i-1)*4+4);
  plot(q_range, Vq/Ry, 'ko-')
  %plot(q_range, -1./(q_range.^2), 'bo-')
  %set(gca,'XScale','log')
  %set(gca,'YScale','log')
  end
  xline(sqrt(3),'r--')
  xline(sqrt(8),'r--')
  xline(sqrt(11),'r--')
  yline(0,'r--')
  xlabel('q (2\pi/a_0)')
  ylabel('V (Ry)')
  box on

%%%%%%%%%%%
% Friedel
  bohr = bohrRadius;
  q_range = linspace(x(1),x(5),101) * 2*pi/a0 * bohrRadius;
  factor = 1;
  a1 = 106.0686 * hartree * (bohrRadius/a0)^2;
  a2 = 2.2278;
  a3 = 0.6060;
  a4 = -1.9720;
  a5 = 5.0;
  a6 = 0.3;
  V_Friedel = a1 * (q_range.^2 - a2)./(exp( a3*(q_range.^2 - a4)) + 1);
  V_Friedel = V_Friedel .* 0.5 .* ( 1+ tanh( (a5-q_range.^2)/a6));
  plot(q_range *a0/bohrRadius *1/(2*pi), V_Friedel/Ry, 'b-')

