% Init concentration
function C = init_C(nb_lat,nb_lat2,nb_lon,distrib)

  C = zeros(nb_lat,nb_lon);
  %if (distrib == 1) % Zonal advection 
  %  for j = 1:20
  %    C(1:90,j) = 80.+j;
  %  end
  %elseif (distrib == 2) % Meridional advection 
  %  C(85:89,:) = 0.5;
  %else 

  % Repartition gaussienne
  sigma2 = 70*70;
  for j = 1:nb_lon
    if (j > 180)
      j2 = j - 360.5;
    else
      j2 = j - 0.5;
    end
    C(1:nb_lat2,j) = exp(-j2*j2*0.5/sigma2);
  end

 % Hemisphere sud
 i_prev = 90;
  for i = nb_lat2+1:nb_lat
    C(i,:) = C(i_prev,:);
    i_prev = i_prev - 1;
  end
end
