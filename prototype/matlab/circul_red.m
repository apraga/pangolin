% Circulation on reduced grid
% circul = 1 : meridional advection
% circul = 2 : zonal advection
% circul = 3 : bidimensionnal advection (hourdin)
% circul = 4 : solid rotation
% circul = 5 : standard test (cv)

function [u,v,dt] = circul_red(t,v_middle,nb_lat,nb_lat2,winds,dlat,period)

  % Winds for solid rotation solide at (lat,lon)
  function result = find_wind_rotation(U_1,is_U,lat,coslat,sinlat,lon,cos_beta0,sin_beta0,cos_alpha,sin_alpha, i,j)
    coslon = cos(lon);
    sinlon = sin(lon);
    % Cartesian
    x = coslat*coslon;
    y = coslat*sinlon;
    z = sinlat;

    % Switching to the new coordinates x1,y1,z1
    x1 = cos_beta0*cos_alpha*x + cos_beta0*sin_alpha*y - sin_beta0*z;
    y1 = -sin_alpha*x + cos_alpha*y;
    z1 = sin_beta0*cos_alpha*x + sin_beta0*sin_alpha*y + cos_beta0*z;

    % Com%pute the spherical coordinates angles
    lat1 = asin(z1); % Dans [-pi/2;pi/2]

    %  Ensure acos domain validity
    frac = x1/cos(lat1);
    prec = 1e-13;
    if (frac > 1. - prec) frac = 1.; end
    if (frac < -1. + prec) frac = -1.; end
    % Special case : longitude in [0,2pi]
    lon1 = acos(frac);

    sinlon1 = y1/cos(lat1);

    if (sinlon1 < 0)
      lon1 = 2*pi - lon1;
    end

    % Speeds computation
    u1 = U_1*cos(lat1);
    % Projection on x1,y1,z1
    u1x = -sin(lon1)*u1;
    u1y = cos(lon1)*u1;

    % Switch to xyz
    ux = cos_beta0*cos_alpha*u1x - sin_alpha*u1y;
    uy = cos_beta0*sin_alpha*u1x + cos_alpha*u1y;
    uz = -sin_beta0*u1x;

    % Projection on (u,v)
    if (is_U)
      result = -ux*sinlon + uy*coslon;
    else
      % v is towards the south so the sign of v is changed
      result = ux*sinlat*coslon + uy*sinlon*sinlat - uz*coslat;
%      if (abs(lat - 89*pi/180.) < 1e-7 && abs(lon - 220*pi/180.) < 1e-7)
%        fprintf('lat lon %.16f %.16f \n', lat, lon);
      %end
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Winds for solid rotation around an axis
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function [u,v] = solid_rotation(u,v,alpha,beta0,nb_lat,v_middle,dlat)
    U_0 = 0.0001;
    radius = 1.;
    % Conversion en degree/s
    U_1 = U_0*180./(pi*radius);
    V_1 = U_0*180./(pi*radius);

    % Position de l'axe de rotation
    % Rotation autour de z, puis rotation autour de y' (notation des angles
    % d'euler)
    cos_beta0 = cos(beta0);
    sin_beta0 = sin(beta0);
    cos_alpha = cos(alpha);
    sin_alpha = sin(alpha);
    coef = pi/180;
    for i=1:nb_lat
      % v est sur les latitudes entieres (1 a 179)
      % u est sur les demi latitudes
      lat_v = (90 - i*dlat)*coef;
      lat_u = (90 - (i-0.5)*dlat)*coef;
      coslat_u = cos(lat_u);
      coslat_v = cos(lat_v);
      sinlat_u = sin(lat_u);
      sinlat_v = sin(lat_v);
      nb_cells = 3*get_nb_cells(i,nb_lat,nb_lat2);
      dlon = 360./nb_cells;
      for j=1:nb_cells
        lon_u = (j-1)*dlon*coef;

%        if (i == 1 && j == 1)
%          fprintf('lat lon for %.15f %.15f', lat_u, lon_u);
%        end 
        u(i,j) = find_wind_rotation(U_1,1,lat_u,coslat_u,sinlat_u,lon_u,cos_beta0,sin_beta0,cos_alpha,sin_alpha, i,j);

      end
      nb_nodes = get_nb_nodes(i,nb_lat,nb_lat2);
      for nb_sec=0:2
        for j=1:nb_nodes
          lon_v = (v_middle(i,j) + nb_sec)*2*pi/3;
          j2 = j + nb_sec*nb_nodes;
          v(i,j2) = find_wind_rotation(U_1,0,lat_v,coslat_v,sinlat_v,lon_v,cos_beta0,sin_beta0,cos_alpha,sin_alpha, i,j);

        end
      end
    end
    % Vitesse nulle au pole sud
    v(nb_lat,:) = 0;
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Bidimensionnelle (hourdin)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function [u,v] = bidim_hourdin(u,v,v_middle,nb_lat,nb_lat2,dlat)
    % Vitesse max des vents (m/s)
    U_0 = 0.0001;
    radius = 1.;
    % Conversion en degree/s
    % ! Les vents sont déjà projetés sur la sphere
    U_1 = U_0*180./(pi*radius);
    V_1 = U_0*180./(pi*radius);
    coef = pi/180. ;

    for i=1:nb_lat
      % v est sur les latitudes entieres (1 a 179)
      % u est sur les demi latitudes
      lat_v = (90 - i*dlat)*coef;
      lat_u = (90 - (i-0.5)*dlat)*coef;
      coslat_u = cos(lat_u);
      coslat_v = cos(lat_v);
      sinlat_u = sin(lat_u);
      sinlat_v = sin(lat_v);
      nb_cells = 3*get_nb_cells(i,nb_lat,nb_lat2);
      dlon = 360./nb_cells;
      for j=1:nb_cells
        lon_u = (j-1)*dlon*coef;
        coslon2_u = cos(lon_u*0.5);
        u(i,j) = U_1*2*coslat_u*sinlat_u*coslon2_u*coslon2_u;
      end
      nb_nodes = get_nb_nodes(i,nb_lat,nb_lat2);
      for nb_sec=0:2
        for j=1:nb_nodes
          lon_v = (v_middle(i,j) + nb_sec)*2*pi/3;
          j2 = j + nb_sec*nb_nodes;
          v(i,j2) = -V_1*coslat_v*cos(lon_v*0.5)*sin(lon_v*0.5);
        end
      end
    end
    % Speed is zero at the south pole
    v(nb_lat,:) = 0;
    % Advection is towards the south so we change the sign of v
    v = -v;
  end

  function [u,v] = meridional(u,v,nb_lat,nb_lat2)
    % Degree/s
    U_0 = 1./160;
    v(:,:) = U_0;
  end

  function [u,v] = zonal(u,v,nb_lat,nb_lat2,dlat)
    % Degree/s
    U_0 = 120./(2*nb_lat2-1)*0.5;
    coef = pi/180;
    for i=1:nb_lat
      u(i,:) = U_0*cos((90 - (i-0.5)*dlat)*coef);
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  nb_cells = 3*get_nb_cells(nb_lat2,nb_lat,nb_lat2);
  nb_nodes = 3*get_nb_nodes(nb_lat2-1,nb_lat,nb_lat2);
  u = zeros(nb_lat,nb_cells);
  v = zeros(nb_lat,nb_nodes);

  if(winds == 1)
    [u,v] = zonal(u,v,nb_lat,nb_lat2,dlat);
    dt = 0.9;

  elseif(winds == 2)
    [u,v] = meridional(u,v,nb_lat,nb_lat2);
    dt = 50;

  elseif(winds == 3)
    [u,v] = bidim_hourdin(u,v,v_middle,nb_lat,nb_lat2,dlat);
    dt = 94.4;
    %dt = 170;

  elseif(winds == 4)
    beta0 = pi*0.5;
    alpha = 0;
    [u,v] = solid_rotation(u,v,alpha,beta0,nb_lat,v_middle,dlat);
    %dt = 43;
    dt = 80;
  elseif(winds == 5)
    %dt = 0.2*period/120.; % 20 lat
    %dt = 0.2*period/120.; % 40 lat
    dt = 0.1*period/120.;
    [u,v] = winds_standard_dv(u,v,t,dt,period,nb_lat,nb_lat2,v_middle,dlat);
  else
    error('other winds not implemented');
  end

end
