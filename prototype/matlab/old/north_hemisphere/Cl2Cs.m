%C_liste vers C_secteur
function [C_s] = Cl2Cs(C_l,nlat,nsec);
  %
  %
  n2=3*nlat*nlat;
  for n=1:n2;
    % n parcourt les mailles
    % i est la latitude
    i=fix(sqrt((n-0.1)/3))+1;
    % j est la longitude de la maille de cette latitude
    j=n-3*(i-1)*(i-1);
    % ns est le numero de secteur
    ns=fix((j-0.1)/(2*i-1))+1;
    % j2 est la longitude
    j2=j-(ns-1)*(2*i-1);
    C(i,j2,ns)=C_l(n);
  end;
 % Pour chaque latitude on recupere toutes les concentrations
  for i=1:nlat;
    for j=1:2*i-1;
      C_s(i,j)=C(i,j,nsec);
    end;
  end;
