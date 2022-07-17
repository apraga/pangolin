%
% Interpollation d'une grille nlat x nlon+1 vers une grille secteur
%
function [u_s] = reg2sec_u(u,mco,dco,nlat,nlon)
%
%
nj1=nlon+1;
for i=1:nlat
    for j=1:nj1;
        X(j)=u(i,j);
    end;
    for nsec=1:3;
    for j=1:2*i-1;
    coo=fix((mco(i,j)-dco(i)*0.5)*120);
    coo=120*(nsec-1)+max(1,coo);
    j1=j+(nsec-1)*(2*i-1);
    u_s(i,j1)=X(coo);
    end;
    end;
end;