%C_secteur vers C_liste
% On stocke d'abord a latitude constante
function [C_l,xlat,xlon] = Cs2Cl(C_s,nlat,nsec);
  % xlat et xlon sont les milieux des mailles 
  % n^2 mailles * 3 secteurs
  n2=3*nlat*nlat;

  C_l = zeros(1,n2);
  xlat = zeros(1,n2);
  xlon = zeros(1,n2);
  %
  dlat=90/nlat;
  for i=1:nlat;
    for j=1:2*i-1;
      dlon=120/(2*i-1);
      % tailles de toutes les mailles stockees jusqu'a la latitude i
      n=3*(i-1)*(i-1);
      np=(nsec-1)*(2*i-1)+j;
      C_l(n+np)=C_s(i,j);
      xlat(n+np)=90-dlat/2-dlat*(i-1);
      xlon(n+np)=120*(nsec-1)+dlon/2+dlon*(j-1);
    end;
  end;
