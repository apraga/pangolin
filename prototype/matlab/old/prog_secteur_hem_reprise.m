% Programme principal
% Advection d'un traceur sur une hemisphere
% avec maillage par secteurs
% suite de la simulation
clear all;
nlat=90;
nlon=360;
dt=20;
nsec=3;
%
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
%
[xlat,xlon] = xlatlon(nlat,nsec);
% Lecture du fichier contenant la distribution du traceur
X=lire_hem(nlat);
C_l=X';
%
C_sto=C_l;
[m0] = masse(C_l,nlat)
%
kk=50;
for k=1:50;
    nit=kk+1
    % Advection 
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
    % Advection v
    %
%    [mass] = masse(C_l,nlat)
    C_l2=C_l*0.;
    for nsec=1:3;
%    [v_s] = reg2sec_v(v,mco,dco,mcoo,dcoo,nlat,nlon,nsec);
    [v_s] = v_ss2v_s(v_ss,nlat,nsec);
    [mC_s] = Cl2Cs(mC_l,nlat,nsec);
    [mC_un_s] = Cl2Cs(mC_un,nlat,nsec);
    [duC_s] = Cl2Cs(duC_l,nlat,nsec);
    [C_s,f_v,C_un] = adv_v7(mC_s,duC_s,mC_un_s,v_s,mco,dco,mcoo,dcoo,nindp,nindm,nlat,nlon,nsec,dt);
    [C_lx,xxlat,xxlon] = Cs2Cl(C_s,nlat,nsec);
    C_l2=C_l2+C_lx;
    end;
    C_l=C_l2;
    kk=kk+1;    
end;
[mass1] = masse(C_l,nlat)
ecr_hem(C_l,nlat);
