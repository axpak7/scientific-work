function y = dydx(y_prev, x)

   % y = x + sqrt(1 + 2 * x^2);
   y = (y_prev + x)/(y_prev - x);
   
end