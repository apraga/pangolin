%      subroutine Circul
function [psi,u,v] = circul (nlat,nlon)
%
pi=2*asin(1.);
dlat=90/nlat;
dlon=360/nlon;
%
% psi: fonction de courant
%
for i=1:nlat+1;
    for j=1:nlon+1;
        psi(i,j)=0.;
    end
end
% initialisation de psi
%
for i=1:nlat;
%    phase(i)=180*pi*(i-1)/(nlat*180)-pi/2;
%    phase(i)=sqrt(phase(i)*phase(i));
    phase(i)=0.;
%    phase(i)=0.
    for j=1:nlon+1;
        psi(i,j)= sin((i-1)*pi/nlat)*cos(2*(j-1)*pi/(nlon)-2*phase(i));
    end
end
%
% Calcul de u et v
%
for i=1:nlat;
    for j=1:nlon+1;
        u(i,j)= (psi(i+1,j)-psi(i,j))/dlat;
    end
end
for i=1:nlat+1;
    coslat=cos((90-(i-1)*dlat)*pi/180);
    for j=1:nlon;
        v(i,j)= -(psi(i,j+1)-psi(i,j))/(dlon*coslat);
    end
end
%v=-v;
%
%u=u*0.;%+0.025;
%v=v*0.;%+0.025;
