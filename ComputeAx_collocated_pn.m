
% ComputeAx
function R=ComputeAx_collocated_pn(X, rr, tt ,delr, delth, r_n, t_n,t,Pneumann)
% int ii r direction,jj theta direction; double P_rr, P_tt, P_r
%X=reshape(X1,r_n,t_n);

R=zeros(r_n,t_n);

for ii=2:r_n-1 %r direction
    for jj=1:t_n % theta direction
        % for theta derivative
        if(jj==1)
            
            P_tt = (X(ii,jj+1)+X(ii,t_n-1)-2*X(ii,jj))/(delth^2);
            
        elseif (jj==t_n)
            P_tt = (X(ii,2)+X(ii,jj-1)-2*X(ii,jj))/(delth^2);
        else
            P_tt = (X(ii,jj+1)+X(ii,jj-1)-2*X(ii,jj))/(delth^2);
        end
        
        % for first r derivatives
        if (ii==2)
            
            P_r=(X(ii+1,jj)-X(ii,jj) + delr*Pneumann(ii,jj) )/(2*delr);% neumann
            P_rr=(X(ii+1,jj)-2*X(ii,jj) + X(ii,jj) -delr*Pneumann(ii,jj))/(delr^2);
            
        elseif (ii==r_n-1)
            % At the farfield
            %if (tt(ii,jj)<pi/6 || tt(ii,jj)>11*pi/6 )
            %   P_r=(0-X(ii-1,jj))/(2*delr);
            % downstream
            %else
            
            P_r=(X(ii,jj)+ delr*Pneumann(ii,jj)-X(ii-1,jj))/(2*delr);% neumann everywhere
            P_rr=(X(ii,jj)+ delr*Pneumann(ii,jj) -2*X(ii,jj)+X(ii-1,jj))/(delr^2);
            %end
            
        else
            P_r=(X(ii+1,jj)-X(ii-1,jj))/(2*delr); %
            P_rr=(X(ii+1,jj)-2*X(ii,jj)+X(ii-1,jj))/(delr^2);
        end
        
        % for second r derivatives
        
        
        R(ii,jj)= P_rr +  (1/(rr(ii,jj)^2))*P_tt + (1/rr(ii,jj))*P_r ;
        
    end
end
%R(r_n-1,1:2)=X(r_n-1,1:2);
%%% R(r_n-3,t_n-3)=X(r_n-3,t_n-3);

%R=R1(:);

end
