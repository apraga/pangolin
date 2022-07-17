prologues := 2; % Necessary for good BB !
filenametemplate "%j-%c.mps"; % for v0.9 or inferior
%outputtemplate := "%j-%c.mps"; % for v1.0 or superior
u:=1cm;
nb_lat2 := 90;
nb_lat := 2*nb_lat2;

% Random
randomseed := 20;

def get_nb_cells(expr i)=
  if (i > nb_lat2) :
    nb := 2*(nb_lat-i+1)-1;
  else:
    nb := 2*i-1;
  fi;
enddef;


% Draw the reduced grid in (lat,lon) coordinates
% Draws only lines
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

def draw_grid_fill(expr n)=
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
      p := (lon,lat)--(lon_next,lat)--(lon_next,lat_next)--(lon,lat_next)--cycle;
      %draw p;
      col := uniformdeviate(1);
      fill p withcolor col*white;
    endfor;
  endfor;
enddef;

% Read comma separated file
vardef read_file(expr s)=
  numeric r[];
  numeric i;
  i = 0;
  forever:
    string line;
    line := readfrom s;
    exitif line = EOF;
    string token, desc;
    k := 0;
    for token = scantokens line:
      if (k = 2):
        r[i] := token;
      fi;
      k := k+1;
    endfor;
    i := i+1;
  endfor;
  r
enddef;

beginfig(1);
  numeric c[];
%  draw_grid_fill(90);
  read_file("list_partitions.txt");
endfig;
end
