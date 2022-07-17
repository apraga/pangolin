%% Analytical solution 
function [u,v] = circul_analytical(nlat,nlon,is_horizontal)
  dlat=90/nlat;
  dlon=360/nlon;

  normV = 0;
  for i = 1:nlat+1
    for j = 1:nlon+1
      lat = (90 - (i-1)/dlat)*pi/180;
      lon = (j-1)/dlon*pi/180;
      if (is_horizontal == 1)
        u(i,j) = 1;
        %u(i,j) = cos(lat); 
        v(i,j) = 0; 
      else
        u(i,j) = 0;
        v(i,j) = 1; 
      end

      %normV = normV + u(i,j)*u(i,j)+v(i,j)*v(i,j);
    end
  end
  normV = (nlat+1)*(nlon+1);
  % Renormalization
  sqrtnorm = sqrt(normV);
  if (is_horizontal == 1)
    u = u / sqrtnorm;
  else
    v = v / sqrtnorm;
  end
end
