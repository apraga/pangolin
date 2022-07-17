%C_secteur vers C_liste
function [C_l,xlat,xlon] = Cs2Cl(C_s,nlat,nsec);
%
%
n2=3*nlat*nlat;
for n=1:n2;
    C_l(n)=0.;
    xlat(n)=0.;
    xlon(n)=0.;
end;
%
dlat=90/nlat;
for i=1:nlat;
    for j=1:2*i-1;
        dlon=120/(2*i-1);
        n=3*(i-1)*(i-1);
        np=(nsec-1)*(2*i-1)+j;
        C_l(n+np)=C_s(i,j);
        xlat(n+np)=90-dlat/2-dlat*(i-1);
        xlon(n+np)=120*(nsec-1)+dlon/2+dlon*(j-1);
    end;
end;
