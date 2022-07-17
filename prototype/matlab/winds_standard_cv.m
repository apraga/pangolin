%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Winds for standard test (convergent version)
% Winds are deformable 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u,v] = winds_standard_cv(u,v,t,dt,period,nb_lat,nb_lat2,v_middle,dlat);
  r = 1.;

  % Earth-like, for testing
  %r = 6.3172e6;
  %eriod = 12*24*3600;
  %dt = 40*60;
  %fprintf('winds between %f %f\n', t, t+dt);

  % m/s -> degree/s
  convert = 360./(2*pi*r);

  % For testing
  %convert = 1;

  U_1 = (10*r/period)*convert;
  U_2 = (2*pi*r/period)*convert;
 
  coef = pi/180;
  const1 = 2*pi*t/period; 
  %int_cos_t = period/(dt*pi)*(sin(pi*(t+dt)/period) - sin(pi*t/period)) ;
  % Half the time step
  cos_t = cos(pi*(t+0.5*dt)/period) ;
  %int_cos_t = cos(pi*t/period);
  for i=1:nb_lat
    lat_v = (90 - i*dlat)*coef;
    lat_u = (90 - (i-0.5)*dlat)*coef;
    coslat_u = cos(lat_u);
    coslat_v = cos(lat_v);
    sin2lat_u = sin(2*lat_u);
    nb_cells = 3*get_nb_cells(i,nb_lat,nb_lat2);
    dlon = 360./nb_cells;
    %if (i == 1)
    %  fprintf('U1 U2 %.10f %.10f \n', U_1, U_2);
    %end
    for j=1:nb_cells
      lon_u = (j-1)*dlon*coef - const1;
      u(i,j)= U_1*sin(lon_u)*sin(lon_u)*sin2lat_u*cos_t + U_2*coslat_u;

    end
    nb_nodes = get_nb_nodes(i,nb_lat,nb_lat2);
    for nb_sec=0:2
      for j=1:nb_nodes
        lon_v = (v_middle(i,j) + nb_sec)*2*pi/3 - const1;
        j2 = j + nb_sec*nb_nodes;
        v(i,j2) = U_1*sin(2*lon_v)*coslat_v*cos_t;

      end
    end
  end
  v = -v;

end


