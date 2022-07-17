
function [] = visu_latlon(C_l,xlat,xlon,nlat);
%
pi=2*acos(0.);
nmax=3*nlat*nlat;
%r=(90-xlat)/90;
for n=1:nmax;
    clat(n)=xlat(n);
    clon(n)=xlon(n);
end;
%
F=TriScatteredInterp(clat',clon',C_l');
%
dmx=0:1:90;
dmy=0:1:360;
[qx,qy]=meshgrid(dmx,dmy);
qz=F(qx,qy);
contourf(qz');