% Advection est-ouest par secteur

function [q_s,dq_s,mCon_s]  =  adv_u_alone(C_s,u_s,nb_lat,nb_lat2,dt)

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Calcul la variation de concentration
   function [delta_f, delta_f1] = variation_C(i,j,q,dq,dlon,dt)
     if (j > 1)
       am = dq(j-1);
       bm = q(j-1)-am*0.5*dlon;
     else % Attention au bord
       am = dq(n);
       bm = q(n)-am*0.5*dlon;
     end

     ap = dq(j);
     bp = q(j)-ap*0.5*dlon;
     
     u_val = u_s(i,j);
     % Flux entre (i,j-1) et (i,j), associé à j-1
     % Si u < 0 : flux de (i,j) vers (i,j-1) (donc negatif)
     if (u_val < 0.)
       dx = -u_val*dt;
       delta_f = -ap*dx*dx*0.5-bp*dx;
       delta_f1 = -dx;
     % Si u > 0 : flux de (i,j-1) vers (i,j) (donc positif)
     else
       dx = u_val*dt;
       delta_f = am*0.5*(2*dlon-dx)*dx+bm*dx; 
       delta_f1 = dx;
     end
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   % Calcul des flux
   dlat = 90/nb_lat2;
   mass = 0;
 
   n = 3*get_nb_mesh(nb_lat2,nb_lat,nb_lat2);
   mC = zeros(1,n);
   mC_un = zeros(1,n);
   q = zeros(1,n);
   dq = zeros(1,n);
   mCon_s = zeros(nb_lat,n); 

   for i = 1:nb_lat;
     nb_mesh = get_nb_mesh(i,nb_lat,nb_lat2);
     n = 3*nb_mesh; % Pour les 3 secteurs

     coslat = cos((90-dlat/2-(i-1)*dlat)*pi/180);
     q(:) = C_s(i,:);
     %
     % Calcul des dérivées
     dlon = 360*coslat/n;
     for j = 2:n-1;
       % Limitation de la pente (van Leer)
       dq(j) = (q(j+1)-q(j-1))*0.5/dlon;
       dC_min = 2*min(abs(q(j+1)-q(j)),abs(q(j-1)-q(j)))/dlon;
       ypm = min(abs(dq(j)),dC_min);
       dq(j) = sign(dq(j))*ypm;
     end;
     % Conditions aux limites périodiques
     dq(1) = (q(2)-q(n))*0.5/dlon;
     dC_min = 2*min(abs(q(2)-q(1)),abs(q(n)-q(1)))/dlon;
     ypm = min(abs(dq(1)),dC_min);
     dq(1) = sign(dq(1))*ypm;
     
     dq(n) = (q(1)-q(n-1))*0.5/dlon;
     dC_min = 2*min(abs(q(1)-q(n)),abs(q(n-1)-q(n)))/dlon;
     ypm = min(abs(dq(n)),dC_min);
     dq(n) = sign(dq(n))*ypm;
      
     for j = 1:n
       [delta_f, delta_f1] = variation_C(i,j,q,dq,dlon,dt);
       f_u(i,j) = delta_f;
       f_u1(i,j) = delta_f1;
     end;

     for j = 1:n-1
       % On additionne les flux aux 2 interfaces j-1/2 et j+1/2
       dx(j) = (-f_u(i,j+1)+f_u(i,j))*dlat;
       dxun(j) = (-f_u1(i,j+1)+f_u1(i,j))*dlat;
       % Mise à jour des masses 
       mCon(j) = dx(j)+q(j)*dlat*dlon;
       mxun(j) = dxun(j)+dlat*dlon;
       %
     end;
     % Traitement de l'autre bord
     dx(n) = (-f_u(i,1)+f_u(i,n))*dlat;
     dxun(n) = (-f_u1(i,1)+f_u1(i,n))*dlat;
     mCon(n) = dx(n)+q(n)*dlat*dlon;
     mxun(n) = dxun(n)+dlat*dlon;
     
     % On reorganise la masse dans les cellules
     for j = 1:n
       q(j) = mCon(j)/mxun(j);
     end;
     
     for j = 2:n-1;
       dq(j) = (q(j+1)-q(j-1))*0.5/dlon;
     end;
     dq(1) = (q(2)-q(n))*0.5/dlon;
     dq(n) = (q(1)-q(n-1))*0.5/dlon;
     
     q_s(i,:) = q(:); 
     mCon_s(i,1:n) = mCon(1:n); 
     dq_s(i,:) = dq(:); 
  end
end
