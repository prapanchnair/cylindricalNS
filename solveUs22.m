function solveUs22(Us)

global Us2;
global rr;
global tt;
global delr;
global delth;
global r_n;
global t_n;
global Re;
global dt;
global k;
b =zeros(size(Us2));


% Define the RHS
for i=2:r_n-1
    for j=1:t_n
        
        if i==r_n-1

              b(i,j) = -Us(i,j)-(dt/Re)*( (1/delr^2)*Us(i+1,j) + (3/rr(i,j))*((1/(2*delr))*Us(i+1,j)));

        else
           b(i,j) = -Us(i,j) ;  
        end
    end
end


  % Biconjugate Gradient Stabilized method %%%%%%%%%%%%%%%%%%%
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%---------------Gonna try matlab bicgstab------------------   
   % Vectors for BICGSTAB
     v0=zeros(size(Us2));
     d0=zeros(size(Us2));
     dd0=zeros(size(Us2));
     x0=eye(size(Us2));%zeros(size(XX));
     %x0(10,10)=1;
     % Scalars
     p0=v0;
     alpha =1;
     w0=1;
     rho0=1;
     
     iter = 0;
     
     Ax0=ComputeAx_collocatedUs2(x0,rr,tt, delr, delth, r_n, t_n,dt,Re);
 
      for i=2:r_n-1
          for j=1:t_n
              dd0(i,j)=b(i,j)-Ax0(i,j);
              d0(i,j)=dd0(i,j);
          end
      end
 %dd0=b-Ax0;
 
 %d0=dd0;
 
     
     Tol=1e-15;
     Res=100;
     
     while Res>Tol,
         rho_i=innerpUV(dd0,d0,r_n,t_n);
         if (rho_i==0) break, end
         if (iter>0)
         beta = (rho_i/rho0)*(alpha/w0);
         p_i = d0 + beta.*(p0 - w0.*v0);
        
         else
             p_i=d0;
            
         end
         v_i = ComputeAx_collocatedUs2(p_i,rr,tt,delr,delth,r_n,t_n,dt,Re);
        
         alpha = rho_i/innerpUV(dd0,v_i,r_n,t_n);
         s = d0 - alpha.*v_i;
        
         t=ComputeAx_collocatedUs2(s,rr,tt,delr,delth,r_n,t_n,dt,Re);
        
         w_i=innerpUV(t,s,r_n,t_n)/innerpUV(t,t,r_n,t_n);
         x_i=x0+alpha.*p_i + w_i.*s;
        
         d_i=s - w_i.*t;
         
         %P= x_i;
         %%%%%%%%
         d0=d_i;
         Res = innerpUV((x_i - x0 ), (x_i -x0),r_n,t_n);
         Res = sqrt ( Res / ((r_n-2)*(t_n)));
         x0=x_i;
         p0=p_i;
         w0=w_i;
         v0=v_i;
        rho0=rho_i;
         
         iter=iter+1;
         
         %fprintf('\n Iter is %d Residual is %0.3g ', iter, Res);
     end
     
     for i=2:r_n-1
         for j=1:t_n
            Us2(i,j) = x_i(i,j);
         end
     end
     
     
     fprintf('\n U converged at Iter %d , Residual is %0.3g at time step %d', iter, Res, k);
end