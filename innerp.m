function d=innerp(A,B,rn,tn)
 d=0;
 for i=2:rn-1
     for j=1:tn
         d=d + A(i,j)*B(i,j);
     end 
 end
 
 end