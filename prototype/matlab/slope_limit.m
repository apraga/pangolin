% Limitation pente (van leer)
% dans le cas horizontal on a l'equivalence
% q_1 <->  q_{i-1} 
% q_2 <->  q_i 
% q_3 <->  q_{i+1} 
function delta = slope_limit(q_1,q_2,q_3,dx)
  sign1 = sign(q_3-q_2);
  sign2 = sign(q_3-q_1);
  sign3 = sign(q_2-q_1);
  delta = 0;
  if (sign1 == sign2 && sign2 == sign3) 
    dC_min = 2*min(abs(q_3-q_2),abs(q_1-q_2));
    dq = 0.5*(q_3-q_1);
    delta = min(abs(dq),dC_min);
    delta = sign(dq)*delta/dx;
  end 

end
