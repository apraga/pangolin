%C_secteur vers C_liste
function [xlat,xlon] = xlatlon(nlat,nsec);
%
%
n2=3*nlat*nlat;
for n=1:n2;
    xlat(n)=0.;
    xlon(n)=0.;
end;
xxlat=xlat;
xxlon=xlon;
%
dlat=90/nlat;
for ns=1:nsec;
for i=1:nlat;
    for j=1:2*i-1;
        dlon=120/(2*i-1);
        n=3*(i-1)*(i-1);
        np=(ns-1)*(2*i-1)+j;
        xlat(n+np)=90-dlat/2-dlat*(i-1);
        xlon(n+np)=120*(ns-1)+dlon/2+dlon*(j-1);
    end;
end;
end;
