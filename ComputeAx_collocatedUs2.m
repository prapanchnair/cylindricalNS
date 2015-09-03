
% ComputeAx
function R=ComputeAx_collocatedUs2(X, rr, tt ,delr, delth, r_n, t_n,dt,Re)
% int ii r direction,jj theta direction; double P_rr, P_tt, P_r
%X=reshape(X1,r_n,t_n);

R=zeros(r_n,t_n);
for ii=2:r_n-1 %r direction
    for jj=1:t_n % theta direction
        % for theta derivative
        if(jj==1)
            
            Uss_tt = (X(ii,jj+1)+X(ii,t_n)-2*X(ii,jj))/(delth^2);
            
        elseif (jj==t_n)
            Uss_tt = (X(ii,1)+X(ii,jj-1)-2*X(ii,jj))/(delth^2);
        else
            Uss_tt = (X(ii,jj+1)+X(ii,jj-1)-2*X(ii,jj))/(delth^2);
        end
        
        % for first r derivatives
        if (ii==2)
            Uss_r=(X(ii+1,jj)-0)/(2*delr);% dirichlet
        elseif (ii==r_n-1)
            % At the farfield
            %     if (tt(ii,jj)<pi/2 || tt(ii,jj)>3*pi/2 ) % neumann
            %         Uss_r=(X(ii,jj)-X(ii-1,jj))/(2*delr);
            %         % downstream
            %    else
            
            Uss_r=(0-X(ii-1,jj))/(2*delr);% Dirichlet
            %     end
            
        else
            Uss_r=(X(ii+1,jj)-X(ii-1,jj))/(2*delr); %
        end
        
        % for second r derivatives
        if (ii==2)
            Uss_rr=(X(ii+1,jj)-2*X(ii,jj) + 0)/(delr^2);
            
        elseif (ii==r_n-1)
            
            %     if (tt(ii,jj)<pi/2 || tt(ii,jj)>3*pi/2 )
            %        Uss_rr=(-1*X(ii,jj)+X(ii-1,jj))/(delr^2);
            
            %   else
            
            
            Uss_rr=(0-2*X(ii,jj)+X(ii-1,jj))/(delr^2);% Dirichlet
            %    end
            
        else
            Uss_rr=(X(ii+1,jj)-2*X(ii,jj)+X(ii-1,jj))/(delr^2);
        end
        
        
        R(ii,jj)= (dt/Re)*( Uss_rr +  (1/(rr(ii,jj)^2))*Uss_tt + (3/rr(ii,jj))*Uss_r ) ...
            +((dt/(Re*rr(ii,jj)^2)) -1)*X(ii,jj);
        
    end
end
%R(r_n-1,1:2)=X(r_n-1,1:2);
%R=R1(:);

end
