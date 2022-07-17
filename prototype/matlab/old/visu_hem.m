
function [] = visu_hem(C_l,xlat,xlon,nlat);
%
pi=2*acos(0.);
nmax=3*nlat*nlat;
r=(90-xlat)/90;
for n=1:nmax;
    clat(n)=r(n)*cos(xlon(n)*pi/180);
    clon(n)=r(n)*sin(xlon(n)*pi/180);
end;
%
F=TriScatteredInterp(clat',clon',C_l');
%
dmx=-1:0.01:1;
dmy=-1:0.01:1;
[qx,qy]=meshgrid(dmx,dmy);
qz=F(qx,qy);
contourf(qz);