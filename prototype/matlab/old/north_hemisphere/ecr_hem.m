function [] = ecr_hem(C_l,nlat);
  %
  fid= fopen ('FI_C_l.txt','w');
  %count = fprintf(fid,'volume mixing ratio\n');
  np=3*nlat*nlat;
  for n=1:np;
    count= fprintf (fid,'  %+10.5e',C_l(n));
  end;

