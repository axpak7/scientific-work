function y = dydx(y_prev, x)

   y = single((y_prev + x)/(y_prev - x));
   
end