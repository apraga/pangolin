function print_maxmin(u,v,s)
  u_max = max(max(u));
  u_min = min(min(u));
  v_max = max(max(v));
  v_min = min(min(v));

  fprintf(strcat(s,' max u %f, max v %f \n'),u_max,v_max);
  fprintf(strcat(s,' min u %f, min v %f \n'),u_min,v_min);
  dx = 2*pi/537;
  dy = pi/180;
  %dt_estim = 0.5/(u_max/dx+v_max/dy);
  dt_estim = 0.5/(abs(u_max)/dx+abs(v_max)/dy);
  fprintf('estimation dt %f\n',dt_estim);

end
