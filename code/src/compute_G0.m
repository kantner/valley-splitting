function G0 = compute_G0(eps, par)
% shift vector for selection rule computation
%  G0 = 2*[ -eps(1,3); -eps(2,3); 1 - eps(3,3)] * 2*pi/par.a0;

  Id = eye(3);
  b1 = 2*pi/par.a0 * [-1,1,1]';
  b2 = 2*pi/par.a0 * [1,-1,1]';

  G0 = (Id - eps)*(b1+b2);

end