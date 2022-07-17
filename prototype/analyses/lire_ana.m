function [lat,lon,lev,A,B,psol,U,V,W,T] = lire_ana(XX)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%ncid= netcdf.open('FMGLOB22+2012012100.nc','NC_NOWRITE');
ncid= netcdf.open(XX,'NC_NOWRITE');
vari=netcdf.inqVarID(ncid,'lat');
lat=netcdf.getVar(ncid,vari);
vari=netcdf.inqVarID(ncid,'lon');
lon=netcdf.getVar(ncid,vari);
vari=netcdf.inqVarID(ncid,'lev');
lev=netcdf.getVar(ncid,vari);
vari=netcdf.inqVarID(ncid,'A_array');
A=netcdf.getVar(ncid,vari);
vari=netcdf.inqVarID(ncid,'B_array');
B=netcdf.getVar(ncid,vari);
vari=netcdf.inqVarID(ncid,'SURFPRESSION');
psol=netcdf.getVar(ncid,vari);
vari=netcdf.inqVarID(ncid,'VENT_ZONAL');
U=netcdf.getVar(ncid,vari);
vari=netcdf.inqVarID(ncid,'VENT_MERIDIE');
V=netcdf.getVar(ncid,vari);
vari=netcdf.inqVarID(ncid,'TEMPERATURE');
T=netcdf.getVar(ncid,vari);
vari=netcdf.inqVarID(ncid,'VITESSE_VERT');
W=netcdf.getVar(ncid,vari);
%end;

