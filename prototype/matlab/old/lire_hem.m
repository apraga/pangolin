function [X] = lire_hem(nlat);
%
fid= fopen ('FI_C_l.txt','r');
%A = fscanf(fid,'%g',[1]);
np=3*nlat*nlat;
%for n=1:np;
       X= fscanf (fid,'%g',[np]);
%end;
