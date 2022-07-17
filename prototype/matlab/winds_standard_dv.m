%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Winds for standard test (divergent version)
% Winds are deformable 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u,v] = winds_standard_dv(u,v,t,dt,period,nb_lat,nb_lat2,v_middle,dlat);
  r = 1.;
  % m/s -> degree/s
  convert = 360./(2*pi*r);
  U_1 = (5*r/period)*convert;
  V_1 = (5*r/(2*period))*convert;
  U_2 = (2*pi*r/period)*convert;
 
  coef = pi/180;
  const1 = 2*pi*t/period; 
  int_cos_t = period/(dt*pi)*(sin(pi*(t+dt)/period) - sin(pi*t/period)) ;
  %int_cos_t = cos(pi*t/period);
  for i=1:nb_lat
    lat_v = (90 - i*dlat)*coef;
    lat_u = (90 - (i-0.5)*dlat)*coef;
    coslat_u = cos(lat_u);
    sinlat_2u = sin(2*lat_u);
    coslat_v = cos(lat_v);
    sin2lat_u = sin(2*lat_u);
    nb_cells = 3*get_nb_cells(i,nb_lat,nb_lat2);
    dlon = 360./nb_cells;
    for j=1:nb_cells
      lon_u = (j-1)*dlon*coef - const1;
      sin_lon_u2 = sin(lon_u*0.5);
      u(i,j) = -U_1*sin_lon_u2*sin_lon_u2*sinlat_2u*coslat_u*coslat_u*int_cos_t + U_2*coslat_u;
    end
    nb_nodes = get_nb_nodes(i,nb_lat,nb_lat2);
    for nb_sec=0:2
      for j=1:nb_nodes
        lon_v = (v_middle(i,j) + nb_sec)*2*pi/3 - const1;
        j2 = j + nb_sec*nb_nodes;
        v(i,j2) = U_1*sin(lon_v)*coslat_v*coslat_v*coslat_v*int_cos_t;
      end
    end
  end
  v = -v;

end


