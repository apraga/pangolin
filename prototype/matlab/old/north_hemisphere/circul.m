%      subroutine Circul
function [u,v] = circul (nlat,nlon)
  %
  pi=2*asin(1.);
  dlat=90/nlat;
  dlon=360/nlon;

 % psi: fonction de courant
 % initialisation de psi

  for i=1:nlat+1;
    for j=1:nlon+1;
      psi(i,j)=0.;
    end
  end
  for i=1:nlat;
    % Non centré au pôle
   phase(i)=0;%pi/4;
    for j=1:nlon+1;
      psi(i,j)= sin((i-1)*pi/nlat)*cos(2*(j-1)*pi/(nlon)-2*phase(i));
    end
  end

  u = zeros(nlat+1,nlon+1);
  v = zeros(nlat+1,nlon+1);
  for i=1:nlat
    u(i,:) = cos((90-i)*pi/180);
  end
  u = u/180.;
  %u = u/norm(u);

  %% Calcul de u et v a partir du potentiel : V = rot psi
  %% u = -1/a d psi/ d lat
  %for i=1:nlat;
  %  for j=1:nlon+1;
  %    u(i,j)= (psi(i+1,j)-psi(i,j))/dlat;
  %  end
  %end
  %%% u = 1/(a cos lat) d psi/ d long
  %for i=1:nlat+1;
  %  coslat=cos((90-(i-1)*dlat)*pi/180);
  %  for j=1:nlon;
  %    v(i,j)= -(psi(i,j+1)-psi(i,j))/(dlon*coslat);
  %  end
  %end

end
