%-------------------------------------------------------------------------------
% Plotting functions on the sphere
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Draw the sphere in orthogonal projection
%-------------------------------------------------------------------------------
% Dot product for vectors
def dot_product(expr a,b)=
  (Xp(a)*Xp(b)+ Yp(a)*Yp(b)+ Zp(a)*Zp(b))
enddef;

% Vectorial product for vectors
def vect_product(expr a,b)=
  (Yp(a)*Zp(b) - Zp(a)*Yp(b), - Xp(a)*Zp(b) + Zp(a)*Xp(b), + Xp(a)*Yp(b) - Yp(a)*Xp(b))
enddef;

% Screen projection
def project(expr V,a,b)=
  (dot_product(V,a),dot_product(V,b))
enddef;

%-------------------------------------------------------------------------------
% Coordinates of a point in spherical coordinates
% phi is the colatitude 
%-------------------------------------------------------------------------------
def coordinates(expr r,theta,phi)=
  (r*cosd(theta)*sind(phi),r*sind(theta)*sind(phi),r*cosd(phi))
enddef;


%-------------------------------------------------------------------------------
% Create (V1,V2,V3) basis and plot it
%-------------------------------------------------------------------------------
def create_basis(expr thetaO,phiO)=
  vector V[];
  V1 = (cosd(thetaO),sind(thetaO),0);
  V2 = (sind(thetaO)*sind(phiO),-cosd(thetaO)*sind(phiO),cosd(phiO));
  V3 = vect_product(V1,V2);

  % Plot contour and vectors
  r = 5*u;
  path cont;
  cont := fullcircle scaled 2r; % Contour
  draw cont;
  % Basis
  z0 = (0,0);
  z1 = project((r,0,0),V1,V2);
  z2 = project((0,r,0),V1,V2);
  z3 = project((0,0,r),V1,V2);
  drawarrow z0--z1;
  drawarrow z0--z2;
  drawarrow z0--z3;
enddef;

%-------------------------------------------------------------------------------
% Search minimal and maximal visible longitude at given colatitude
% set nb_intersect and lon_min,lon_max on the stack
% lat_max and lat_min are set only for 2 intersections
%-------------------------------------------------------------------------------
def get_visible_lon(expr r,phi)=
  % The hidden part is found as the intersection between the plane orthogonal to
  % V3 and the parallel : % V3x*x + V3y*y + V3z*z = 0
  % with
  % X= r*cos(theta) sin(phi), 
  % Y = r*sin(theta) sin(phi), 
  % Z = r*cos(phi)

  alpha := Xp(V3);
  beta := Yp(V3);
  gamma := Zp(V3);
  Z := r*cosd(phi);
  a := (alpha**2+beta**2);
  b := 2*alpha*gamma*Z;
  c := (gamma*Z)**2 - (r*beta*sind(phi))**2;
  % X is root of a 2nd order polynomial
  if (a = 0): 
    X1 := -c/b;
    Y1 := -alpha*X1-gamma*Z;
    nb_intersect := 1;
  else:
    delta := b*b - 4*a*c;
    if (delta < 0):
      nb_intersect := 0;
    elseif (delta = 0):
      X1 := -b/(2*a);
      % We deduce Y from X
      Y1 := (-alpha*X1-gamma*Z)/beta;
      nb_intersect := 1;
    else:
      X1 := (-b - sqrt(delta))/(2*a);
      X2 := (-b + sqrt(delta))/(2*a);
      Y1 := (-alpha*X1-gamma*Z)/beta;
      Y2 := (-alpha*X2-gamma*Z)/beta;
      nb_intersect := 2;
    fi;
  fi;
  lat_min := 0;
  lat_max := 0;
  if (nb_intersect = 2):
    lat_min := angle(X1,Y1);
    lat_max := angle(X2,Y2);
    if (lat_min > lat_max):
      tmp := lat_max;
      lat_max := lat_min;
      lat_min := tmp;
    fi;
  fi;
enddef;

%-------------------------------------------------------------------------------
% Draw a parallel at colatitude phi
% Hidden part is drawn with dashed line
%-------------------------------------------------------------------------------
def draw_parallel(expr r,phi)=
  front := 1; %1 if we are in front of (x,y,z)
  get_visible_lon(r,phi);
  if ((thetaO > 0) and (thetaO < 180)):
    front := -1;
  fi;
  %% Draw part of the parallel in solid line and the hidden part in dashed line
  %% (or invisible according to view_hidden)

  if (nb_intersect = 2):
    sign := 1;
    if (lat_min < 0):
      sign := -1;
    fi;
    % sign*front < 0 if it is opposed to the camera
    draw_arc_colat(lat_min,lat_max,phi,sign*front);
    draw_arc_colat(lat_max,360+lat_min,phi,-sign*front);
  % Draw circles near the pole
  else:
    sign := 1;
    if (phi > 90):
      sign := -1;
    fi;
    draw_arc_colat(0,360,phi,sign*front);
  fi;

enddef;

%-------------------------------------------------------------------------------
% Draw an arc at given colatitude
%-------------------------------------------------------------------------------
def draw_arc_colat(expr t_init,t_last,colat,is_visible)=
  path p;
  dt := 10;
  p = project(coordinates(r,t_init,colat),V1,V2)
  for t=t_init+dt step dt until t_last:
    .. project(coordinates(r,t,colat),V1,V2)
  endfor
    .. project(coordinates(r,t_last,colat),V1,V2);
  if (is_visible > 0):
    draw p withcolor blue;
  elseif (view_hidden > 0):
    draw p dashed evenly withcolor red;
  fi;
enddef;

%-------------------------------------------------------------------------------
% Search minimal and maximal visible colatitude at given longitude
% set nb_intersect and lat_min,lat_max on the stack
% lat_max and lat_min are set only for 2 intersections
%-------------------------------------------------------------------------------
def get_visible_lat(expr lon)=
  % Compute the visible part
  alpha := Xp(V3)*cosd(lon);
  beta := Yp(V3)*sind(lon);
  gamma := Zp(V3);
  lat_min := angle((alpha+beta),-gamma);
  lat_max := lat_min + 180;
  if (lat_min < 0):
    tmp := lat_max;
    lat_max := 360 + lat_min;
    lat_min := tmp;
  fi;
enddef;

%-------------------------------------------------------------------------------
% Draw an arc at given longitude
% Draws only the visible part
%-------------------------------------------------------------------------------
def draw_arc_lon(expr t_init,t_last,lon,is_visible)=
  path p;
  dt := 10;
  get_visible_lat(lon);
 
  if ((t_init < lat_max) and (t_last > lat_min)):
    p = project(coordinates(r,lon,lat_min),V1,V2)
    for t=lat_min+dt step dt until lat_max-dt:
      .. project(coordinates(r,lon,t),V1,V2)
    endfor
      .. project(coordinates(r,lon,lat_max),V1,V2);
  fi;
  %if (is_visible > 0):
  draw p withcolor blue;
  %elseif (view_hidden > 0):
  %  draw p dashed evenly withcolor red;
  %fi;
enddef;

%-------------------------------------------------------------------------------
% Number of cells at latitude i
%-------------------------------------------------------------------------------
def get_nb_cells(expr i)=
  if (i > nb_lat2) :
    nb := 2*(nb_lat-i+1)-1;
  else:
    nb := 2*i-1;
  fi;
enddef;

%-------------------------------------------------------------------------------
% Draw the reduced grid in (lat,lon) coordinates
% Draws only lines
%-------------------------------------------------------------------------------
def draw_grid(expr n)=
  width = 10*u;
  height = 10*u;
  step_y = height/n;
  %show n;
  for i=0 upto n-1:
    cur_i := height - i*step_y;
    cur_i_next := cur_i - step_y;
    draw (0,cur_i)--(width,cur_i);
    nb := get_nb_cells(i+1);
    step_x := width/nb;
    for j=0 upto nb:
      cur_j := j*step_x;
      draw (cur_j,cur_i)--(cur_j,cur_i_next);
    endfor;
  endfor;
  draw (0,0)--(width,0);
enddef;

%-------------------------------------------------------------------------------
% Draw partition with colors stored in file in (lat,lon)
%-------------------------------------------------------------------------------
def draw_partition_fill(expr n)=
  width = 10*u;
  height = 10*u; 
  dlat = height/n;
  path p;
 
  for i=0 upto n-1:
    lat := height - i*dlat;
    lat_next := lat - dlat;
    nb := get_nb_cells(i+1);
    dlon := width/nb;
    for j=0 upto nb-1:
      lon := j*dlon;
      lon_next := lon + dlon;

      read_line("list_partitions.txt"); 
      p := (lon,lat)--(lon_next,lat)--(lon_next,lat_next)--(lon,lat_next)--cycle;
      fill p withcolor (r,g,b);
      draw p;
    endfor;
  endfor;
enddef;

%-------------------------------------------------------------------------------
% Read only one line in comma separated file
% and store it into (r,g,b)
% Faster than using large array
%-------------------------------------------------------------------------------
def read_line(expr s)=
  string line;
  string token;
  line := readfrom s;
  numeric k;
  k = 0;
  exitif line = EOF;
  for token = scantokens line:
    % Read only rgb values
    if (k = 2):
      r := token;
    elseif (k = 3):
      g := token;
    elseif (k = 4):
      b := token;
    fi;
    k := k+1;
  endfor;
enddef;


