%      subroutine Circul
% distrib = 1 : u = cst, v = 0
% distrib = 2 : u = 0, v = 1
function [u,v] = circul (nb_lat,nb_lon,radius,distrib)
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
   % v > 0 correspond à des vents vers le pole sud
   U_0 = 0.001;
   % Conversion de m/s en rad/s
   V_1 = U_0/(2*pi*radius);
    for i=1:nb_lat
      lat = (90 - i*dlat)*pi/180;
      sinlat2 = sin(2*lat);
      coslat = cos(lat);
      U_1 = U_0/(2*pi*radius*coslat);
      if (coslat == 0)
        U_1 = 1.;
      end
      for j=1:nb_lon
        lon = j*dlon*pi/180;
        %lon = (j*dlon-180)*pi/180;
        cos2lon = cos(lon*0.5);
        sinlon = sin(lon);
        u(i,j)= U_1*sinlat2*cos2lon*cos2lon;
        v(i,j)= -V_1*0.5*coslat*sinlon;
      end
    end
    %U_0 = norm(u)+norm(v);
    %u = u/U_0;
    %v = v/U_0;
    x = zeros(nb_lat,nb_lon);
    y = zeros(nb_lat,nb_lon);
    for i=1:nb_lat
      x(i,:) = [0:nb_lon-1];
    end
    for j=1:nb_lon
      y(:,j) = [90:-1:-89];
    end
    quiver(x,y,u,v,5);
% L'advection se fait vers le sud donc on inverse v
   v = -v;

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
