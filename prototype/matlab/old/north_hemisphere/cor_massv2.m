%
% Correction de la composante v de la vitesse
% si v constant

function [v_ss] = cor_massv2(v,mcoo,dcoo,indm,nlat,nlon)
  %
  %
  for i=1:nlat;
    for n=1:3;
      resv(i)=0.;
      for j=1:(4*i-1)
        coo=fix((mcoo(i,j)-dcoo(i,j))*120);
        coo=120*(n-1)+max(1,coo);
        v_ss(i,j,n)=v(i,coo);
        %fprintf('i,j = %d %d gives coo %d \n',i,j,coo);
      end;
      for j=1:4*i-1;
        resv(i)=resv(i)+v_ss(i,j,n)*dcoo(i,j)/3;
      end;
      %fprintf('resv %f \n',resv(i));
      resv(i) = resv(i) / (4*i-1);
      tmp = 0.055 / (3*(4*i-1));
      %fprintf('resv %f - tmp %f \n',resv(i),tmp);

      for j=1:4*i-1;
        %if (v_ss(i,j,n) ~= 0)
        %  fprintf('vss %f - resv %f \n',v_ss(i,j,n),resv(i));
        %end
        v_ss(i,j,n)=v_ss(i,j,n)-resv(i);
        %if (v_ss(i,j,n) ~= 0)
        %  fprintf('i,j = %d %d %f\n',i,j,v_ss(i,j,n));
        %end
      end;
    end;
    lol = sum(v_ss(i,:,:))
  end;
