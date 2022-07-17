% Concentration directly on reduced grid
% Types of distribution
% distrib = 1 : zonal 
% distrib = 2 : meridional 
% distrib = 3 : gaussian hills
% distrib = 4 : cosine bells
% distrib = 5 : slotted cylinders
% distrib = 6 : cosine bells correlated
% default : gaussian

function q_s = init_C_red(nb_lat,nb_lat2,distrib,dlat)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  function [x,y,z] = get_cartesian(lat,lon,r)
    x = r*cos(lon)*cos(lat);
    y = r*sin(lon)*cos(lat);
    z = r*sin(lat);
  end

  % Gaussian repartition 
  function [q_s] = gaussian(q_s,nb_lat,nb_lat2)
    sigma2 = 70*70;
    for i = 1:nb_lat
      nb_cells = 3*get_nb_cells(i,nb_lat,nb_lat2);
      dlon = 360./nb_cells;
      for j = 1:nb_cells
        lon = (j-0.5)*dlon;
        if (lon > 180.)
          lon = 360. - lon;
        end
        q_s(i,j) = exp(-lon*lon*0.5/sigma2);
      end
    end
  end

  % Gaussian hills
  function [q_s] = gaussian_hills(q_s,lat1,lat2,lon1,lon2,nb_lat,nb_lat2,dlat)
    hmax = 0.95;
    b = 5;
    r = 1.;
    coef = pi/180;
    [x1,y1,z1] = get_cartesian(lat1,lon1,r);
    [x2,y2,z2] = get_cartesian(lat2,lon2,r);
    for i = 1:nb_lat
      nb_cells = 3*get_nb_cells(i,nb_lat,nb_lat2);
      lat = (90 - (i-0.5)*dlat)*coef;
      dlon = 360./nb_cells;
      for j = 1:nb_cells
        lon = (j-0.5)*dlon*coef;
        [x,y,z] = get_cartesian(lat,lon,r);
        dist1 = (x-x1)*(x-x1) + (y-y1)*(y-y1) + (z-z1)*(z-z1);
        dist2 = (x-x2)*(x-x2) + (y-y2)*(y-y2) + (z-z2)*(z-z2);
        h1 = hmax*exp(-b*dist1);
        h2 = hmax*exp(-b*dist2);
        q_s(i,j) = h1 + h2;
      end
    end
  end

  % Compute the great distance on the sphere between (lon1,lat1) and (lon,lat)
  % Input in radians
  function r = great_distance(lon1,cos_lat1,sin_lat1,lon,cos_lat,sin_lat,r)
    r = r*acos(sin_lat1*sin_lat + cos_lat1*cos_lat*cos(lon-lon1));
  end

  % Cosine bells distribution
  function [q_s] = cosine_bells(is_linear,q_s,lat1,lat2,lon1,lon2,nb_lat,nb_lat2,dlat)
    hmax = 1;
    b = 0.1;
    c = 0.9;
    a_phi = -0.8;
    b_phi = 0.9;
    r = 1.;
    r_base = r*0.5;
    coef = pi/180;
    cos_lat1 = cos(lat1);
    sin_lat1 = sin(lat1);
    cos_lat2 = cos(lat2);
    sin_lat2 = sin(lat2);
    for i = 1:nb_lat
      nb_cells = 3*get_nb_cells(i,nb_lat,nb_lat2);
      lat = 90 - (i-0.5)*dlat;
      dlon = 360./nb_cells;
      cos_lat = cos(lat*coef);
      sin_lat = sin(lat*coef);
      for j = 1:nb_cells
        lon = (j-0.5)*dlon*coef;
        r1 = great_distance(lon1,cos_lat1,sin_lat1,lon,cos_lat,sin_lat,r);
        cur = b;
        if (r1 < r_base)
          h1 = hmax*0.5*(1+cos(pi*r1/r_base));
          cur = b + c*h1;
        else
          r2 = great_distance(lon2,cos_lat2,sin_lat2,lon,cos_lat,sin_lat,r);
          if (r2 < r_base)
            h2 = hmax*0.5*(1+cos(pi*r2/r_base));
            cur = b + c*h2;
          end
        end
        if (is_linear)
          q_s(i,j) = cur;
        else
          q_s(i,j) = a_phi*cur*cur + b_phi;
        end
      end
    end
  end

  function q_s = slotted_cylinders(q_s,lat1,lat2,lon1,lon2,nb_lat,nb_lat2,dlat)
    r = 1.;
    r_base = r*0.25;
    c = 0.9;
    b = 0.1;
    coef = pi/180;
    cos_lat1 = cos(lat1);
    sin_lat1 = sin(lat1);
    cos_lat2 = cos(lat2);
    sin_lat2 = sin(lat2);

    for i = 1:nb_lat
      nb_cells = 3*get_nb_cells(i,nb_lat,nb_lat2);
      lat = (90 - (i-0.5)*dlat)*coef;
      dlon = 360./nb_cells;
      cos_lat = cos(lat);
      sin_lat = sin(lat);
      for j = 1:nb_cells
        lon = (j-0.5)*dlon*coef;
        r1 = great_distance(lon1,cos_lat1,sin_lat1,lon,cos_lat,sin_lat,r);
        r2 = great_distance(lon2,cos_lat2,sin_lat2,lon,cos_lat,sin_lat,r);

        q_s(i,j) = b;
        if (r1 <= r_base)
          if (r*abs(lon - lon1) >= r_base/6)
            q_s(i,j) = c;
          elseif (r*(lat - lat1) < -5*r_base/12)
            q_s(i,j) = c;
          end
        elseif (r2 <= r_base)
          if (r*abs(lon - lon2) >= r_base/6)
            q_s(i,j) = c;
          elseif (r*(lat - lat2) > 5*r_base/12)
            q_s(i,j) = c;
          end
        end
      end
    end
  end

  function zonal_distrib(i1,i2,nb_lat,nb_lat2)
    for i = i1:i2
      nb_cells = 3*get_nb_cells(i,nb_lat,nb_lat2);
      dlon = 360./nb_cells;
      for j = 1:nb_cells
        q_s(i,j) = 1.;
      end
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

  nb_cells_max = 3*get_nb_cells(nb_lat2,nb_lat,nb_lat2);
  q_s = zeros(nb_lat,nb_cells_max);
  lat1 = 0;
  lat2 = 0;
  lon1 = 5*pi/6;
  lon2 = 7*pi/6;

  % q is a ratio 
  % For meridional advection
  if (distrib == 1 ) 
    zonal_distrib(41,50,nb_lat,nb_lat2);

  % For zonal advection
  elseif (distrib == 2) 
    width = 10.;
    for i = 1:nb_lat
      nb_cells = 3*get_nb_cells(i,nb_lat,nb_lat2);
      dlon = 360./nb_cells;
      for j = 1:nb_cells
        lon = (j-0.5)*dlon;
        if (lon < width)
          q_s(i,j) = 1.;
        end
      end
    end

  elseif (distrib == 3) 
    q_s = gaussian_hills(q_s,lat1,lat2,lon1,lon2,nb_lat,nb_lat2,dlat);

  elseif (distrib == 4) 
    q_s = cosine_bells(1,q_s,lat1,lat2,lon1,lon2,nb_lat,nb_lat2,dlat);

  elseif (distrib == 5) 
    q_s = slotted_cylinders(q_s,lat1,lat2,lon1,lon2,nb_lat,nb_lat2,dlat);

  elseif (distrib == 6) 
    q_s = cosine_bells(0,q_s,lat1,lat2,lon1,lon2,nb_lat,nb_lat2,dlat);

  % For solid rotation
  elseif (distrib == 7 ) 
    zonal_distrib(81,100,nb_lat,nb_lat2);

  else
    q_s = gaussian(q_s,nb_lat,nb_lat2);
  end
end
