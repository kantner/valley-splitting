function G0 = compute_G0(par)
% shift vector for selection rule computation
  G0 = 2*[ -par.eps(1,3); -par.eps(2,3); 1 - par.eps(3,3)] * 2*pi/par.a0;

end