%-------------------------------------------------------------------------------
% Contains Metapost routines for drawing the partitions
%-------------------------------------------------------------------------------
prologues := 2; % Necessary for good BB !
filenametemplate "%j-%c.mps"; % for v0.9 or inferior
%outputtemplate := "%j-%c.mps"; % for v1.0 or superior
warningcheck := 0; % No warning about numbers > 4096
u:=1cm;
nb_lat2 := 90;
nb_lat := 2*nb_lat2;

beginfig(1);
  input lib_plot
  let vector=color;
  let Xp = redpart;
  let Yp = greenpart;
  let Zp = bluepart;

  % Theta must be > 180
  thetaO = 210;
  % Phi0 must be < 0
  phiO = -10;
  create_basis(thetaO,phiO);
  view_hidden := -1;

  %for colat=1 step 1 until 180:
  %  draw_parallel(r/u,colat); %Dividing by u avoids overflow
  %endfor;
  for lon=0 step 10 until 360:
    draw_arc_lon(0,180,lon,1);
  endfor;

endfig;

end
