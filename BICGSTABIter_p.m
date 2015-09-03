function [p_int] = BICGSTABIter_p(u,r,Nr,Np,hr,hp,f)

%ip = Nr-1 ; jp = Np/2 ;

STOP = false ;

BICGEPS = 1.0e-16 ;

Temp=ComputeAX_p(u,r,Nr,Np,hr,hp) ;

% Initial vector r_0 = b - Ax_0, and r0* = r_0

for i = 2:(Nr-1)
    for j =1:Np
        
        Uj(i,j) = u(i,j) ;
        rj(i,j) = f(i,j) - Temp(i,j) ;
        r0_star(i,j) = rj(i,j) ;
        
        %{
            
            if (i == ip) && (j == jp)
              
                Uj(i,j) = 0.0 ;
                rj(i,j) = 0.0 ;
                r0_star(i,j) = 0.0 ;
                         
            end
            
        %}
        
    end
end

BICG_ITER = 0 ; norm = 0.0 ;

while ((BICG_ITER < 10000 ) && (~STOP == 1))
    
    % compute rhoj = (r0, r0*)
    
    rhoj = 0.0 ;
    for i = 2:(Nr-1)
        for j =1:Np
            rhoj = rhoj + rj(i,j)*r0_star(i,j) ;
        end
    end
    
    if ( sqrt(rhoj/((Nr-1)*(Np+1))) < BICGEPS )
        
        STOP = true ;
        
    else
        
        if ( BICG_ITER == 0 )
            
            for i = 2:(Nr-1)
                for j = 1:Np
                    pj(i,j) = rj(i,j); % p0 = r0
                end
            end
            
        else
            betaj = (rhoj/rhoj_Minus)*(alphaj/omegaj) ;
            
            for i = 1:(Nr-1)
                for j = 1:Np
                    pj(i,j) = rj(i,j) + betaj*(pj(i,j) - omegaj*Var(i,j));
                end
            end
        end
        
        %Solve for Upstar, Vpstar from Ku* = u...., where K is the preconditioning matrix
        
        for i = 2:(Nr-1)
            for j = 1:Np
                % No preconditioning
                pstar(i,j) = pj(i,j) ;
            end
        end
        
        % compute vj = A*pstar
        
        for i = 2:(Nr-1)
            for j = 1:Np
                
                %{
                            
                     if (i == ip) && (j == jp)
                         
                         u(i,j) = 0.0 ;
                         
                     else
                         
					u(i,j) = pstar(i,j) ;
                    
                     end
                %}
                
                
                u(i,j) = pstar(i,j) ;
                
            end
        end
        
        Temp=ComputeAX_p(u,r,Nr,Np,hr,hp) ;
        
        for i = 2:(Nr-1)
            for j = 1:Np
                Var(i,j) = Temp(i,j) ;
            end
        end
        H1 = 0.0 ;
        for i = 2:(Nr-1)
            for j = 1:Np
                H1 = H1 + Var(i,j)*r0_star(i,j) ;
            end
        end
        alphaj = rhoj/H1 ;
        
        % find sj
        
        for i = 2:(Nr-1)
            for j = 1:Np
                sj(i,j) = rj(i,j) - alphaj*Var(i,j) ;
            end
        end
        
        % Solve for Upstar, Vpstar from Ku* = u...., where K is the preconditioning matrix
        
        for i = 2:(Nr-1)
            for j = 1:Np
                % No preconditioning
                sstar(i,j) = sj(i,j) ;
            end
        end
        norm = 0.0 ;
        for i = 2:(Nr-1)
            for j = 1:Np
                norm = norm + sstar(i,j)*sstar(i,j) ;
            end
        end
        norm = sqrt(norm/((Nr-1)*(Np+1))) ;
        if( norm < BICGEPS)
            STOP = true ; % if ||s||_2 is small x_i = x_{i-1} + alphai*p_i
            for i = 2:(Nr-1)
                for j = 1:Np
                    
                    %{
                                
                                if (i == ip) && (j == jp)
                                    
                                    Uj(i,j) = 0.0 ;
                                    
                                else
                                    
	                		        Uj(i,j) = Uj(i,j)+ alphaj*pstar(i,j) ;
                                    
                                end
                                
                    %}
                    Uj(i,j) = Uj(i,j)+ alphaj*pstar(i,j) ;
                end
            end
        else
            
            % compute t = As
            
            for i = 2:(Nr-1)
                for j = 1:Np
                    
                    %{
                                
                                if (i == ip) && (j == jp)
                                    
                                    u(i,j) = 0.0 ;
                                    
                                else
                                    
						            u(i,j) = sstar(i,j) ;
                        
                                end
                                
                    %}
                    u(i,j) = sstar(i,j) ;
                end
            end
            
            Temp=ComputeAX_p(u,r,Nr,Np,hr,hp) ;
            
            H1 = 0.0 ; H2 = 0.0 ;
            for i = 2:(Nr-1)
                for j = 1:Np
                    H1 = H1 + Temp(i,j)*sj(i,j) ;
                    H2 = H2 + Temp(i,j)*Temp(i,j) ;
                end
            end
            omegaj = H1/H2;
            
            % find xj
            norm = 0.0 ;
            for i = 2:(Nr-1)
                for j = 1:Np
                    H1 = (alphaj*pstar(i,j) + omegaj*sstar(i,j)) ;
                    %{
                        
                        if (i == ip) && (j == jp)
                            
                            Uj(i,j) = 0.0 ;
                            
                        else
                            
                            Uj(i,j) = Uj(i,j) +  H1 ;
                            
                        end
                    %}
                    Uj(i,j) = Uj(i,j) +  H1 ;
                    norm = norm + H1*H1 ;
                end
            end
            norm = sqrt(norm/((Nr-1)*(Np+1))) ;
            
            if(norm < BICGEPS)
                STOP = true ;
            end
            
            % find rjplusone
            
            for i = 2:(Nr-1)
                for j = 1:Np
                    rj(i,j) = sj(i,j) - omegaj*Temp(i,j);
                end
            end
            rhoj_Minus = rhoj ;
        end
    end
    
    BICG_ITER = BICG_ITER + 1 ;
    
end

for i = 2:(Nr-1)
    for j = 1:Np
        p_int(i,j) = Uj(i,j) ;
    end
end

fprintf('\n Pressure converged at Iter %d, ', BICG_ITER);

end

