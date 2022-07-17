% Calcule la distribution des vents

function [u,v] = circul(nb_lat,nb_lon,distrib)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function view_vectorfield(nb_lat,nb_lon,u,v)
    x = [1:nb_lon];
    y = [1:nb_lat];
    y = 90 -y;
    quiver(x,y,u,v,2);
  end

  % Vents pour rotation solide autour d'un axe
  function [u,v] = solid_rotation(alpha,beta0,nb_lat,nb_lon)
    U_0 = 0.0001;
    % Conversion en degree/s
    U_1 = U_0*180/(pi*radius);
    V_1 = U_0*180/(pi*radius);

    % Position de l'axe de rotation
    % Rotation autour de z, puis rotation autour de y' (notation des angles
    % d'euler)
    cos_beta0 = cos(beta0);
    sin_beta0 = sin(beta0);
    cos_alpha = cos(alpha);
    sin_alpha = sin(alpha);
    dlat = 1.;
    dlon = 1.;
    k = 1;
    for i=1:nb_lat
      lat = (90.5 - i*dlat)*pi/180;
      for j=1:nb_lon
        lon = (j-0.5)*dlon*pi/180;

        % Coordonnees x,y,z
        x = cos(lat)*cos(lon);
        y = cos(lat)*sin(lon);
        z = sin(lat);

        % Passage dans le nouveau systeme de coordonnees x1,y1,z1
        x1 = cos_beta0*cos_alpha*x + cos_beta0*sin_alpha*y - sin_beta0*z;
        y1 = -sin_alpha*x + cos_alpha*y;
        z1 = sin_beta0*cos_alpha*x + sin_beta0*sin_alpha*y + cos_beta0*z;
        % Calcul des angles
        lat1 = asin(z1);
        % Longitude dans [0,2pi] donc etude a part
        lon1 = acos(x1/cos(lat1));
        % Si sin(lon1) < 0
        sinlon1 = y1/cos(lat1);
        if (sinlon1 < eps)
          lon1 = 2*pi - lon1;
        end
        %fprintf('lat avant apres %f %f \n',lat,lat1);
        % Calcul des vitesses
        u1 = U_1*cos(lat1);
        % Projection sur x1,y1,z1
        u1x = -sin(lon1)*u1;
        u1y = cos(lon1)*u1;

        % Passage dans xyz
        ux = cos_beta0*cos_alpha*u1x - sin_alpha*u1y;
        uy = cos_beta0*sin_alpha*u1x + cos_alpha*u1y;
        uz = -sin_beta0*u1x;

        % Projection sur (u,v)
        u(i,j) = -ux*sin(lon) + uy*cos(lon);
        v(i,j) = -ux*sin(lat)*cos(lon) - uy*sin(lon)*sin(lat) + uz*cos(lat);
      end
    end
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  nb_lat2 = nb_lat/2;
  dlat=180/nb_lat;
  dlon=360/nb_lon;

  u = zeros(nb_lat,nb_lon);
  v = zeros(nb_lat,nb_lon);
  radius = 1.;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Zonal
  if (distrib == 1)
    for i=1:nb_lat
      u(i,:) = cos((90.5-i)*pi/180);
    end
    U_0 = 0.0001;
    % Conversion en degree/s
    U_0 = U_0*180/(pi*radius);
    u = u*U_0;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif(distrib == 2)
    U_0 = 0.0001;
    V_1 = U_0*180/(pi*radius);
    v(:,:) = V_1;
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Bidimensionnelle (hourdin)
  elseif(distrib == 3)
    % Vitesse max des vents (m/s)
    U_0 = 0.0001;
    % Conversion en degree/s
    % ! Les vents sont déjà projetés sur la sphere
    U_1 = U_0*180/(pi*radius);
    V_1 = U_0*180/(pi*radius);

    for i=1:nb_lat
      % v est sur les latitudes entieres (1 a 179)
      % u est sur les demi latitudes
      lat_v = (90 - i*dlat)*pi/180;
      lat_u = (90.5 - i*dlat)*pi/180;
      coslat_u = cos(lat_u);
      coslat_v = cos(lat_v);
      sinlat_u = sin(lat_u);
      sinlat_v = sin(lat_v);
      for j=1:nb_lon
        lon_v = (j-0.5)*dlon*pi/180;
        % Premiere vitesse est en longitude = 0
        lon_u = (j-1)*dlon*pi/180;
        coslon2_u = cos(lon_u*0.5);
        coslon2_v = cos(lon_v*0.5);
        sinlon2_v = sin(lon_v*0.5);
        u(i,j)= U_1*2*coslat_u*sinlat_u*coslon2_u*coslon2_u;
        v(i,j)= -V_1*coslat_v*coslon2_v*sinlon2_v;
      end
    end
    % Vitesse nulle au pole sud
    v(nb_lat,:) = 0;
    % L'advection se fait vers le sud donc on inverse v
    v = -v;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Bidimensionnelle
  elseif(distrib == 4)
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

  else
    beta0 = pi*0.25;
    alpha = 0;
    [u,v] = solid_rotation(alpha,beta0,nb_lat,nb_lon);

    % Vitesse nulle au pole sud
    v(nb_lat,:) = 0;
    %view_vectorfield(nb_lat,nb_lon,u,v);
    %% L'advection se fait vers le sud donc on inverse v
    v = -v;
  end
end
