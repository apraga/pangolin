clear all;
close all;
XX='FMGLOB22+2012012100.nc';
nlat=90;
nlon=180;
nlev=60;
p00=101325.;
kapa=2/7;
[lat,lon,lev,A,B,psol,U,V,W,T] = lire_ana(XX);
for j=1:nlat;
    for i=1:nlon;
        for k=1:nlev;
            P(i,j,k)=A(k)*p00+B(k)*psol(i,j);
        end;
    end;
end;
pmoy=(A+B)*p00/100;
% mid troposphere:
%Teta0=330.;
% tropopause:
%Teta0=370.;
% basse stratosphere:
%Teta0=500.;
% moyenne stratosphere:
Teta0=850.;
% haute stratosphere
%Teta0=1200.;
%
for j=1:nlat;
for i=1:nlon;
%for k=1:nlev;
for k=1:50;    
            TETA(k)=T(i,j,k)*exp(kapa*(log(p00)-log(P(i,j,k))));
            UU(k)=U(i,j,k);
            VV(k)=V(i,j,k);
            TT(k)=T(i,j,k);
            PP(k)=P(i,j,k);
            niveau(k)=k;
end;
%type='spline';
type='linear';
%type='nearest';
UTET (i,j)= interp1(TETA, UU,Teta0,type);
VTET (i,j)= interp1(TETA, VV,Teta0,type);
PTET (i,j)= interp1(TETA, PP,Teta0,type);
TTET (i,j)= interp1(TETA, TT,Teta0,type);
ETET (i,j)= UTET(i,j)*UTET(i,j) + VTET(i,j)*VTET(i,j);
end;
end;
tratet;