
function [] = gnu_write(iter,C_l,xlat,xlon,nlat);
  %
  name = strcat ('../plot/tmp/',num2str(iter),'.dat');
  fid= fopen (name,'w');
  np=3*nlat*nlat;
  for n=1:np;
    x = cos(xlon(n)*pi/180)*cos(xlat(n)*pi/180);
    y = sin(xlon(n)*pi/180)*cos(xlat(n)*pi/180);
    z = sin(xlat(n)*pi/180);
    count= fprintf (fid,' %f %f %f %+10.5e \n \n',x,y,z,C_l(n));
  end;
  fclose(fid);
