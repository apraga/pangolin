% East-west advection

function [q_s,dq_s,mCon_s,mCon_un_s] = adv_u_alone(C_s,u_s,nb_lat,nb_lat2,dt,dlat)

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Variation of concentration
   function [delta_f, delta_f1] = variation_C(i,j,q,dq,dlon,dt,n)
     if (j > 1)
       j_prev = j-1;
     else % Beware of the sides
       j_prev = n;
     end

     am = dq(j_prev);
     bm = q(j_prev)-am*0.5*dlon;

     ap = dq(j);
     bp = q(j)-ap*0.5*dlon;
     
     u_val = u_s(i,j);
%     if (i == 80 && j == 120) 
%       fprintf('zwinds %.16f \n', u_val);
%     end
     % Flux entre (i,j-1) et (i,j), associé à j
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
   
   % Derivatives
   function [dq] = compute_dq(dq,q,n,dlon)
     for j = 2:n-1;
       % Slope limitation (van Leer)
       dq(j) = slope_limit(q(j-1),q(j),q(j+1),dlon);
%       if (i == 80 && j == 120) 
%         fprintf('prev, cur, next %.16f %.16f %.16f\n', q(j-1), q(j), q(j+1));
%         fprintf('dq %.16f\n', dq(j));
%       end

     end;
     % Periodic boundary conditions
     dq(1) = slope_limit(q(n),q(1),q(2),dlon);
     dq(n) = slope_limit(q(n-1),q(n),q(1),dlon);
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   mass = 0;
 
   nb_cells = 3*get_nb_cells(nb_lat2,nb_lat,nb_lat2);
   mC = zeros(1,nb_cells);
   mC_un = zeros(1,nb_cells);
   q = zeros(1,nb_cells);
   dq = zeros(1,nb_cells);
   mCon_s = zeros(nb_lat,nb_cells); 
   flux = zeros(nb_cells,1); 
   flux_un = zeros(nb_cells,1); 

   for i = 1:nb_lat;
     nb_cells = 3*get_nb_cells(i,nb_lat,nb_lat2); % Pour les 3 secteurs

     coslat = cos((90-dlat/2-(i-1)*dlat)*pi/180);
     q(:) = C_s(i,:);

     % Calcul des dérivées
     dlon = 360*coslat/nb_cells;
     dq = compute_dq(dq,q,nb_cells,dlon);
     dS = dlat*dlon;
      
     for j = 1:nb_cells
       [delta_f, delta_f1] = variation_C(i,j,q,dq,dlon,dt,nb_cells);
       % f_{j-1/2} est stocke en j
       flux(j) = delta_f;
       flux_un(j) = delta_f1;
     end;

     for j = 1:nb_cells-1
       % On additionne les flux aux 2 interfaces j et j+1 pour q_{i+1/2}
       dx = (-flux(j+1)+flux(j))*dlat;
       dxun = (-flux_un(j+1)+flux_un(j))*dlat;
       if (i == 80 && j == 119) 
       %  fprintf('flux prev %.16f %.16f\n', flux(j), flux_un(j));
         %`fprintf('dlon %.16f\n', dlon);
         %fprintf('q %.16f\n', q(j));
         %fprintf('slope %.16f\n', dq(j)*dlon);
         %fprintf('qneighb %.16f\n', q(j+1));
         %fprintf('slopeneighb %.16f\n', dq(j+1)*dlon);
         %fprintf('flux next %.16f %.16f\n', flux(j+1), flux_un(j+1));
       end

       % Mise à jour des masses 
       mCon(j) = dx + q(j)*dS;
       mxun(j) = dxun + dS;
     end;
     % Traitement de l'autre bord
     dx = (-flux(1)+flux(nb_cells))*dlat;
     dxun = (-flux_un(1)+flux_un(nb_cells))*dlat;
     mCon(nb_cells) = dx + q(nb_cells)*dS;
     mxun(nb_cells) = dxun + dS;
     
     % On reorganise la masse dans les cellules
     for j = 1:nb_cells
       q(j) = mCon(j)/mxun(j);
%       if (i == 81 && j == 483)
%          fprintf('q after zonal%.16f \n', q(j));
%        end
     end
     
     for j = 2:nb_cells-1;
       dq(j) = (q(j+1)-q(j-1))*0.5/dlon;
     end;
     dq(1) = (q(2)-q(nb_cells))*0.5/dlon;
     dq(nb_cells) = (q(1)-q(nb_cells-1))*0.5/dlon;
%     if (i == 80) 
%         fprintf('nbcells %d\n',nb_cells);
%         fprintf('for slope %.16f %.16f\n',q(nb_cells-1), q(1));
%       end
     
     q_s(i,:) = q(:); 
     mCon_s(i,1:nb_cells) = mCon(1:nb_cells); 
     mCon_un_s(i,1:nb_cells) = mxun(1:nb_cells); 
     dq_s(i,:) = dq(:); 
  end
end
