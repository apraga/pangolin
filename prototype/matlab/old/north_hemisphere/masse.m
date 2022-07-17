%
% advection nord-sud par secteur
%
function [mass] = masse(C,nlat)
%
%
% Calcul des flux
%
dlat=90/nlat;
pi=2*acos(0.);
mass=0;
for i=1:nlat;
    coslat=cos((90-dlat/2-(i-1)*dlat)*pi/180);
    jl=3*(2*i-1);    
    dlon=360*coslat/jl;
    for j=1:jl;
        ind=3*(i-1)*(i-1)+j;
        % concentration = masse / surface 
        mass=mass+C(ind)*dlat*dlon;
        xun(j)=1;
    end;
end;
    %
