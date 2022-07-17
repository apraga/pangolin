% Programme principal
% Advection d'un traceur sur une hemisphere
% avec maillage par secteurs
%
% Definition de la grille
nlat=90;
nlon=360;
% Pas de temps
dt=20;
%
%%%%%%%%%%%%%%%
% Init sur une bande
C = zeros(nlat,nlon);
C(3,:)=1;
%%%%%%%%%%%%%%%

%
% Initialisation par liste et grille secteur
%
[mco,dco] = grille(nlat);
[mcoo,dcoo,isec,nindp,nindm]=sub_grille(nlat);
%
% Initialisation de la circulation
%[u,v] = circul (nlat,nlon);
[u,v] = circul_analytical (nlat,nlon);
[u_s] = reg2sec_u(u,mco,dco,nlat,nlon);
% Correction de u et v pour conservation de la masse
[v_ss] = no_cor_mass_v(v,mcoo,dcoo,nlat,nlon);
u_ss = u_s;
% Pas de conservation de la masse
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

gnu_write(0,C_l,xlat,xlon,nlat);

% Stockage concentration initiale et calcul de la masse
C_sto=C_l;
[m0] = masse(C_l,nlat)
% Integration temporelle 
kk=0;
for k=1:1000
  [m0] = masse(C_l,nlat)
  nit=kk+1;
  % Advection nord-sud de v
  C_l2=C_l*0.;
  for nsec=1:3;    
    [v_s] = v_ss2v_s(v_ss,nlat,nsec);
    [C_s,f_v,C_un] = adv_v_alone(C_l,v_s,mco,dco,mcoo,dcoo,nindp,nindm,nlat,nlon,nsec,dt);
    [C_lx,xxlat,xxlon] = Cs2Cl(C_s,nlat,nsec);
    C_l2=C_l2+C_lx;
  end;
  C_l=C_l2;

  kk=kk+1
  gnu_write(k,C_l,xlat,xlon,nlat);
end;
% Fin de l'elvolution temporelle
% diagnostic masse
[mass1] = masse(C_l,nlat)

