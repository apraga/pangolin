prologues := 2; % Necessary for good BB !
filenametemplate "%plot_grid.mps"; % for v0.9 or inferior
%outputtemplate := "%j-%c.mps"; % for v1.0 or superior
u:=5mm;

% Read comma separated file
% type is an integer for the type of file
def read_file(expr s,type)=
  numeric r[];
  numeric i;
  forever:
    string line;
    line := readfrom s;
    exitif line = EOF;
    string token;
    i := 0;
    for token = scantokens line:
      r[i] := token*u;
      i := i+1;
    endfor;
    if (type = 1):
      plot_cell_partition(r);
    elseif (type = 2):
      plot_neighbours(r,i);
    elseif (type = 3):
      plot_cell(r);
    else:
      show "Type of file not recognized !";
      exitif true;
    fi;
  endfor;
enddef;

% Plot a cell with its value
% Coord contains the cell value, then 4 points (x_i,y_i) defining the corners
% of the cell
def plot_cell_partition(suffix coord)=
  numeric j;
  path p;
  % Input in (lat,lon) so we switch
  p := (coord[2],coord[1]);
  for j = 3 step 2 until i-1:
    p := p--(coord[j+1],coord[j]);
  endfor
  p := p--cycle;
  fill p withcolor colormap[coord[0] mod nb_colors];
  draw p;
enddef;

% Plot a cell with a color corresponding to its partition
% Coord contains the partition id, then 4 points (x_i,y_i) defining the corners
% of the cell
def plot_cell_partition(suffix coord)=
  numeric j;
  path p;
  % Input in (lat,lon) so we switch
  p := (coord[2],coord[1]);
  for j = 3 step 2 until i-1:
    p := p--(coord[j+1],coord[j]);
  endfor
  p := p--cycle;
  fill p withcolor colormap[coord[0] mod nb_colors];
  draw p;
enddef;

% Link a partition to its neighbours
% Format : x0 y0 [x1 y1 [x2 y2]] 
def plot_neighbours(suffix coord,i)=
  numeric  k;
  path p;
  pair origin_i, origin_tmp, cur, vec;
  pickup pencircle scaled 2;
  interim ahangle := 60;
  interim ahlength := 15;
  origin_i := (coord[1],coord[0]);
  draw origin_i withpen pencircle scaled 10;
  for k = 3 step 2 until i:
    cur := (coord[k],coord[k-1]);
    vec := cur - origin_i;
    % Scale a bit for readability
    cur := origin_i + 0.9*vec;
    origin_tmp := origin_i + 0.1*vec;
    p :=origin_tmp--cur ;
    drawarrow p ;
  endfor;
enddef;

% Define a colormap
def write_colormap(expr n)=
  nb_colors = 16;
  colormap[0] = (1.00, 0.00, 0.00); % Red
  colormap[1] = (0.00, 1.00, 0.00); % Green        
  colormap[2] = (1.00, 1.00, 0.00); % Yellow       
  colormap[3] = (0.00, 0.00, 1.00); % Blue         
  colormap[4] = (1.00, 0.00, 1.00); % Magenta      
  colormap[5] = (0.00, 1.00, 1.00); % Cyan         
  colormap[6] = (1.00, 0.50, 0.20); % Orange       
  colormap[7] = (0.30, 0.55, 0.00); % Olive        
  colormap[8] = (0.72, 0.47, 0.47); % Dark pink    
  colormap[9] = (0.33, 0.33, 0.81); % Sea blue     
  colormap[10] = (1.00, 0.63, 0.63); % Pink         
  colormap[11] = (0.62, 0.44, 0.65); % Violet       
  colormap[12] = (0.60, 0.80, 0.70); % Pale green   
  colormap[13] = (0.47, 0.20, 0.00); % Brown        
  colormap[14] = (0.00, 0.68, 0.68); %Turquoise    
  colormap[15] = (0.81, 0.00, 0.40); % Purple     
enddef;
    
beginfig(1);
  color colormap[];
  write_colormap(0);
  read_file("../output/grid.dat",1);
  read_file("../output/neighbours.dat",2);
endfig;
end
