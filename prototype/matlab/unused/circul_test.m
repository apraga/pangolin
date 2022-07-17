%      subroutine Circul
% distrib = 1 : u = cst, v = 0
% distrib = 2 : u = 0, v = 1
function [u,v] = circul (nb_lat,nb_lon,distrib)
  %
  pi=2*asin(1.);
  nb_lat2 = nb_lat/2;
  dlat=180/nb_lat;
  dlon=360/nb_lon;

  % psi: fonction de courant
  % initialisation de psi

  u = zeros(nb_lat,nb_lon);
  v = zeros(nb_lat,nb_lon);
  % Zonal
  if (distrib == 1)
    for i=1:nb_lat
      u(i,:) = cos((90-i)*pi/180);
    end
    % On ne veut pas une vitesse trop grande par rapport a la taille d'une 
    % cellule (min = 2pi/537)
    U_0 = 180.;
    u = u/U_0;

  % Meridional
  elseif(distrib == 2)
    v(:,:) = 1.;
    v = v/norm(v);

  % Bidimensionnelle (hourdin)
  elseif(distrib == 3)
   radius = 1.;%3600000.;
   % Vitesse max des vents (m/s)
   U_0 = 1;
   % Conversion en rad/s
   %U_0 = U_0*360/(2*pi*radius);
   % Normalisation
   V_1 = U_0/(2*radius);

   %phi0 = asin(1./sqrt(3));
   %U_1 = U_0/(radius*pi);%*sin(2*phi0)*cos(phi0));
   U_1 = U_0/(radius*pi);

   % Pas de vitesse aux poles
   % On va de la latitude 1 a 179
   %for i=1:nb_lat
   %  % v est sur les latitudes entieres
   %  % u est sur les demi latitudes
   %  lat_v = (90 - i*dlat)*pi/180;
   %  lat_u = (90.5 - i*dlat)*pi/180;
   %  coslat_u = cos(lat_u);
   %  coslat_v = cos(lat_v);
   %  sinlat_u = sin(lat_u);
   %  sinlat2_u = sin(2*lat_u);
   %  sinlat_v = sin(lat_v);
   %  for j=1:nb_lon
   %    lon_u = (j-0.5)*dlon*pi/180;
   %    lon_v = j*dlon*pi/180;
   %    coslon2_u = cos(lon_u*0.5);
   %    coslon2_v = cos(lon_v*0.5);
   %    sinlon2_u = sin(lon_u*0.5);
   %    sinlon2_v = sin(lon_v*0.5);
   %    sinlon_v = sin(lon_v);
   %    u(i,j)= U_1*2*coslat_u*sinlat_u*coslon2_u*coslon2_u;
   %    v(i,j)= -V_1*sinlat_v*coslon2_v*sinlon2_v;
   %  end
   %end
   for i=1:nb_lat
     lat = (90 - i*dlat)*pi/180;
     coslat = cos(lat);
     sinlat = sin(lat);
     for j=1:nb_lon
       lon = j*dlon*pi/180;
       coslon2 = cos(lon*0.5);
       sinlon2 = sin(lon*0.5);
       u(i,j)= U_1*2*coslat*sinlat*coslon2*coslon2;
       v(i,j)= -V_1*coslat*coslon2*sinlon2;
     end
   end
   
   x = [1:nb_lon];
   y = [1:nb_lat];
   y = 90 -y;
   quiver(x,y,u,v,2);
   % L'advection se fait vers le sud donc on inverse v
   v = -v;
   u = u*pi/180;
   v = v*pi/180;


  % Bidimensionnelle
  else
    psi = zeros(nb_lat+1,nb_lon+1);
    for i=1:nb_lat2
      % Non centré au pôle
      phase(i)=0;%pi/4;
      for j=1:nb_lon+1;
        psi(i,j)= sin((i-1)*pi/nb_lat2)*cos(2*(j-1)*pi/(nb_lon)-2*phase(i));
      end
    end

    % Calcul de u et v a partir du potentiel : V = rot psi
    %u = -1/a d psi/ d lat
    for i=1:nb_lat;
      for j=1:nb_lon;
        u(i,j)= (psi(i+1,j)-psi(i,j))/dlat;
      end
    end
    % u = 1/(a cos lat) d psi/ d long
    for i=1:nb_lat+1;
      coslat=cos((90-(i-1)*dlat)*pi/180);
      for j=1:nb_lon;
        v(i,j)= -(psi(i,j+1)-psi(i,j))/(dlon*coslat);
      end
    end
  end

end
