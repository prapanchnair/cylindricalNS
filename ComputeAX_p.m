function [Temp] = ComputeAX_p(u,r,Nr,Np,hr,hp)

%ip = Nr-1 ; jp = Np/2 ;

for i = 2:(Nr-1)
		for j =1:Np 
            
			if i==2
                
				u_r = (u(i+1,j) - u(i,j))/(2.0*hr) ;
				u_rr = (u(i+1,j) - u(i,j))/(hr*hr) ;
                
            elseif i==(Nr-1)
                
				u_r = (u(i,j) - u(i-1,j))/(2.0*hr) ;
				u_rr = (-u(i,j) + u(i-1,j))/(hr*hr) ;
                
            else
                
                u_r = (u(i+1,j) - u(i-1,j))/(2.0*hr) ;
				u_rr = (u(i+1,j) - 2.0*u(i,j) + u(i-1,j))/(hr*hr) ;
                
            end
			
			if j==1 
                
                u_pp = (u(i,j+1) - 2.0*u(i,j) + u(i,Np))/(hp*hp) ; 
                
            elseif j==Np 
                
                u_pp = (u(i,1) - 2.0*u(i,j) + u(i,j-1))/(hp*hp) ; 
                
            else
                
                u_pp = (u(i,j+1) - 2.0*u(i,j) + u(i,j-1))/(hp*hp) ;
                
            end
            

			%{
            if i==ip && j==jp
                
                Temp(i,j) = 0 ;
                
            else
                
                Temp(i,j) = u_rr + (1.0/r(i))*u_r + (1.0/(r(i)*r(i)))*u_pp ;
                
            end
            
            %}
            
            Temp(i,j) = u_rr + (1.0/r(i))*u_r + (1.0/(r(i)*r(i)))*u_pp ;
            
        end
end

end

