% Simple copie

function [v_ss] = no_cor_mass_v(v,mcoo,dcoo,nlat,nlon)
  %
  %
  for i=1:nlat;
    for j=1:nlon;
      X(j)=v(i+1,j);
    end;
    resv(i)=0.;
    for n=1:3;
      for j=1:4*i-1;
        coo=fix((mcoo(i,j)-dcoo(i,j)*0.5)*120);
        coo=120*(n-1)+max(1,coo);
        j1=j;%+(nsec-1)*(4*i-1);
        v_ss(i,j1,n)=X(coo);
      end;
    end;
  end
end
