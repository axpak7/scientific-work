function y = dydx(x, y_prev)

   y = single((y_prev + x)/(y_prev - x));
   
end