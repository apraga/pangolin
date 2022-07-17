function [X] = gnu_read(iter,nlat);
  name = strcat ('../plot/tmp/',num2str(iter),'.dat');
  fid= fopen (name,'r');
  X= fscanf (fid,'%*g %*g %*g %g ');
  fclose(fid);
