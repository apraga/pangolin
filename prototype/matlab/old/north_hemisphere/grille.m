%
% calcul des coordonnees des secteurs
% com contient le numero des mailles
% dco contient le pas entre les mailles
function [com,dco] = grille(nlat)
  %
  %
  nsec=2*nlat-1;
  for n=1:nlat;
    dco(n)=1./(2*n-1);
    for i=1:2*n-1;
      com(n,i)=(i-0.5)*dco(n);
    end;
  end;
  %
