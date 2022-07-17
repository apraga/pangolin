% Init concentration
function C = init_C(nb_lat,nb_lat2,nb_lon,distrib)

  C = zeros(nb_lat,nb_lon);
  if (distrib == 1) % Zonal advection 
    for j = 1:20
      C(1:90,j) = 80.+j;
    end
  elseif (distrib == 2) % Meridional advection 
    C(85:89,:) = 0.5;
  else 
    for j = 1:nb_lon/2+1
      C(1:nb_lat2,j) = j;
    end
    for j = nb_lon/2+2:nb_lon
      C(1:nb_lat2,j) = nb_lon - j+2;
    end
  end

 % Hemisphere sud
 i_prev = 90;
  for i = nb_lat2+1:nb_lat
    C(i,:) = C(i_prev,:);
    i_prev = i_prev - 1;
  end
end
