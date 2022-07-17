% Programme principal
% Advection d'un traceur sur une hemisphere
% avec maillage par secteurs
%
clear all;
%
% Definition de la grille
nlat=90;
nlon=360;
% Pas de temps
dt=20;
%
% Initialisation du traceur sur une grille latitude longitude
% reguliere
for i=1:nlat,
  for j=1:nlon;
    C(i,j)=i;
  end;
end;
%
% Initialisation par liste et grille secteur
%
[mco,dco] = grille(nlat);
[mcoo,dcoo,isec,nindp,nindm]=sub_grille(nlat);
%
% Initialisation de la circulation
[psi,u,v] = circul (nlat,nlon);
[u_s] = reg2sec_u(u,mco,dco,nlat,nlon);
% Correction de u et v pour conservation de la masse
[v_ss] = cor_massv(v,mcoo,dcoo,nlat,nlon);
[u_ss,um,div,div2m] = cor_massu(u_s,v_ss,dcoo,nindm,nindp,nlat);
% Initialisation du traceur dans chaque secteur
nsec=1;
[C_s,dC_s]= reg2sec(C,mco,dco,nlat,nlon,nsec);
[C_l,xlat,xlon] = Cs2Cl(C_s,nlat,nsec);

%
for nsec=2:3;
  [C_s,dC_s]= reg2sec(C,mco,dco,nlat,nlon,nsec);
  [C_lx,xxlat,xxlon] = Cs2Cl(C_s,nlat,nsec);
  C_l=C_l+C_lx;
  xlat=xlat+xxlat;
  xlon=xlon+xxlon;
end;
% Stokage concentration initiale et calcul de la masse
C_sto=C_l;
[m0] = masse(C_l,nlat)
% Integration temporelle 
kk=0;
for k=1:1500
  nit=kk+1
  % Advection est-ouest
  [mC_l,mC_un,duC_l] = adv_u5(C_l,u_ss,mco,dco,nlat,nlon,dt);
  % Lissage pole
  %   mCpole=(mC_l(1)+mC_l(2)+mC_l(3))/3;
  %   mC_l(1)=mCpole;
  %   mC_l(2)=mCpole;
  %   mC_l(3)=mCpole;
  %   mC_unp=(mC_un(1)+mC_un(2)+mC_un(3))/3;
  %   mC_un(1)=mC_unp;
  %   mC_un(2)=mC_unp;
  %   mC_un(3)=mC_unp;
  %
  % Advection v nord-sud par secteur
  %
  C_l2=C_l*0.;
  for nsec=1:3;    
    [v_s] = v_ss2v_s(v_ss,nlat,nsec);
    [mC_s] = Cl2Cs(mC_l,nlat,nsec);
    [mC_un_s] = Cl2Cs(mC_un,nlat,nsec);
    [duC_s] = Cl2Cs(duC_l,nlat,nsec);
    [C_s,f_v,C_un] = adv_v7(mC_s,duC_s,mC_un_s,v_s,mco,dco,mcoo,dcoo,nindp,nindm,nlat,nlon,nsec,dt);
    %C_s(1:3,1:5)
    [C_lx,xxlat,xxlon] = Cs2Cl(C_s,nlat,nsec);
    C_l2=C_l2+C_lx;
  end;
  C_l=C_l2;
  %C_l(1:3*3*3)
  if (mod(k,50) == 0)
    ncl_write2(k,C_l,nlat,nlon)
  end
  kk=kk+1;    
end;
% Fin de l'elvolution temporelle
% diagnostic masse
[mass1] = masse(C_l,nlat)
% Ecriture du la concentration finale
%ecr_hem(C_l,nlat);
