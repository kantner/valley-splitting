function [y] = log_cosh(x)

  % regular
    y = log(cosh(x));

  % |x| large
    x_high = 710;
    idx    = abs(x) > x_high;
    y(idx) = abs(x(idx)) - log(2);
end