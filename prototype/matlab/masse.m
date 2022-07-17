% Compute the total mass and check for divergences
function [mass] = masse(C_s,nb_lat,nb_lat2,m0,check_type,to_display,dlat)

  mass = 0;
  err = 0.1;
  for i=1:nb_lat;
    coslat = cos((90-dlat/2-(i-1)*dlat)*pi/180);
    nb_cells = 3*get_nb_cells(i,nb_lat,nb_lat2);
    dlon = 360*coslat/nb_cells;
    dS = dlon*dlat;
    for j=1:nb_cells;
      % concentration = masse / surface 
      mass = mass+C_s(i,j)*dS;
    end
  end
  if (to_display == 1)
    fprintf('mass = %.16f, diff=%.16e \n',mass, mass-m0);
  end

  if (check_type == 1)
    if (mass < 0)
      error('divergence après advection zonale !');
    end
  elseif (check_type == 2)
    if (abs(mass - m0) > err)
      error('divergence de la masse !');
    end
  end
end
