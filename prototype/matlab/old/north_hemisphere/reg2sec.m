%
% Interpollation d'une grille nlat x nlon vers une grille secteur
%
function [C_s,dC_s] = reg2sec(C,mco,dco,nlat,nlon,nsec)
  %
  %
  dlat=90/nlat;
  nj1=nlon/3;
  nj2=(nsec-1)*nj1+1;
  for i=1:nlat;
    % Selectionne un seul secteur
    for j=1:nj1;
      X(j)=C(i,nj2+j-1);
    end;
    nppsec=nj1/(2*i-1);
    % On garde seulement les valeurs au milieu de la maille
    for j=1:2*i-1;
      coo=fix(mco(i,j)*120+0.5);
      coo=max(1,coo);
      C_s(i,j)=X(coo);
    end;
  end;
  % Calcul des derivées (2e ordre)
  % d C-s / d latitude
  %
  dC_s(1,1)=0.;
  for i=2:nlat-1;
    for j=1:2*i-1;
      % Numero de la maille donne directement !
      coo1=fix(mco(i,j)/dco(i-1)+0.5);
      coo1=max(1,coo1);
      coo2=fix(mco(i,j)/dco(i+1)+0.5);
      coo2=max(1,coo2);
      dC_s(i,j)=(C_s(i+1,coo2)-C_s(i-1,coo1))*0.5/dlat;    
    end;
  end;
  % Attention au bord
  for j=1:2*nlat-1;
    coo1=fix(mco(nlat,j)/dco(nlat-1)+0.5);
    coo1=max(1,coo1);
    coo2=fix(mco(nlat,j)/dco(nlat)+0.5);
    coo2=max(1,coo2);
    dC_s(nlat,j)=(C_s(nlat,coo2)-C_s(nlat-1,coo1))/dlat;    
  end;

