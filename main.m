% Numerical Solution to Navier Stokes Equation in 2D Polar Coordinates
% Second Order accuracy throughout
% Uses Biconjugate Gradient Stabilized scheme for Matrix Solution.
% Author: Prapanj G R
% Department of Mechanical Engineering, Indian Institute of Science,
% Bangalore, India.

% r, t coordinates. r- radial, t - theta.
% t varies from 0 rad (downstream 0 angle) to 2pi radians. Anticlockwise is
% positive

function main()
clear all; clc; format long;

global Us2;
global Vs2;
global rr;
global tt;
global delr;
global delth;
global r_n;
global t_n;
global Re;
global dt;
global k;

delr    =0.075;
delth   =pi/75; %32
r       =[1:delr:20];
t       =[0:delth:(2*pi-delth)];
r_size  = size(r); r_n= r_size(2);
t_size  = size(t); t_n= t_size(2);
[rr,tt] = ndgrid(r,t);

for i=1:r_n
    for j=1:t_n
        XX(i,j)=rr(i,j)*cos(tt(i,j));
        YY(i,j)=rr(i,j)*sin(tt(i,j));
    end
end
for i=1:r_n
    for j=1:t_n+1
        if j==t_n+1
            gridX(i,j) = XX(i,1);
            gridY(i,j) = YY(i,1);
        else
            gridX(i,j)=XX(i,j);
            gridY(i,j)=YY(i,j);
        end
    end
end

fid = fopen('grid.xyz','w');
fprintf(fid, '%d \t %d\n', r_n, t_n+1);
fprintf(fid,'%f\n',gridX);
fprintf(fid,'%f\n',gridY);
fclose(fid);
writenum= 50;
% In grid, the first argument is in r direction
% 2nd argument is in theta direction.

% let Ux =  10 at the upstream. Then U = -Ux* cos(pi - t)
% V = -Ux*sin(pi - t)

% Initial U and V
U = zeros(size(XX));
V = zeros(size(XX));
Ux = zeros(size(XX));
Uy = zeros(size(XX));
Us = zeros(size(XX));
Vs = zeros(size(XX));
Us2 = zeros(size(XX));
Vs2 = zeros(size(XX));
G1 = zeros(size(XX));
G2 = zeros(size(XX));
deldotU = zeros(size(XX));
P = zeros(size(XX));


UxBc=1.0;
for j=1:t_n
    for i=2:r_n
        
        U(i,j)=UxBc*(1 - (1/r(i))^2)*cos(tt(i,j));
        V(i,j)=-UxBc*(1 + (1/r(i))^2)*sin(tt(i,j));
       
    end
end

% Reynolds number
Re=150;
dt = 0.005;

%  Time iterations
for k=1:40000
    
    % Convection
    for i=2:r_n-1
        for j=1:t_n
            
            
            convU1 = (1/(2*delr))*( U(i+1,j)^2 - U(i-1,j)^2 );
            convV1 = (1/(2*delr))*(U(i+1,j)*V(i+1,j) - U(i-1,j)*V(i-1,j) );
            convV3 = 2*U(i,j)*V(i,j)/rr(i,j);
            convU3 =  (U(i,j)^2 - V(i,j)^2)/rr(i,j); %%%%%% edited
                        
            if (j == 1)  %Theta boundaries
                
                convU2 = (1/rr(i,j))*(1/(2*delth))*(U(i,j+1)*V(i,j+1) - U(i,t_n)*V(i,t_n));
                convV2 = (1/rr(i,j))*(1/(2*delth))*(V(i,j+1)^2 - V(i,t_n)^2);
                            
            elseif (j==t_n)
                
                convU2 = (1/rr(i,j))*(1/(2*delth))*(U(i,1)*V(i,1) - U(i,j-1)*V(i,j-1));
                convV2 = (1/rr(i,j))*(1/(2*delth))*(V(i,1)^2 - V(i,j-1)^2);
                
            else
                
                convU2 = (1/rr(i,j))*(1/(2*delth))*(U(i,j+1)*V(i,j+1) - U(i,j-1)*V(i,j-1));
                convV2 = (1/rr(i,j))*(1/(2*delth))*(V(i,j+1)^2 - V(i,j-1)^2);
            end
            
           Us(i,j) = U(i,j) - dt*( convU1+convU2+convU3); 
           Vs(i,j) = V(i,j) - dt*( convV1+convV2 + convV3); 
            
        end
    end
     
    sum=0;
    for j=1:t_n
        
        sum= sum+U(r_n,j);
    end
    fprintf('\n Sum after convection = %f ',sum);
    
    for j=1:t_n

            Us(r_n,j) = U(r_n,j);
            Vs(r_n,j) = V(r_n,j);
            Us2(r_n,j) = U(r_n,j);
            Vs2(r_n,j) = V(r_n,j);
        
    end
    
    solveUs22(Us);
    solveVs22(Vs);
    
    for i=2:r_n-1
        for j=1:t_n
            
            U(i,j)=Us2(i,j);
            V(i,j)=Vs2(i,j);
            
        end
    end
    
    sum=0;
    for j=1:t_n
        
        sum= sum+U(r_n,j);
    end
    fprintf('\n Sum after diffusion = %f ',sum);

 %   U(1,:)=0;V(1,:) =0;
 
    for i = 2: r_n-1
        for j = 1:t_n
            if (i==2) % wall
                term1= (U(i+1,j) - 0)/(2*delr) ;
                %        term1= (1/(4*delr))*( (rr(i+1,j)+rr(i,j))*(U(i+1,j)+U(i,j)) - (rr(i-1,j) +rr(i,j))*(0 +U(i,j))  );
            elseif (i==r_n-1)
                
                term1= (U(i+1,j) - U(i-1,j))/(2*delr) ;
                
            else
                
                term1= (U(i+1,j) - U(i-1,j))/(2*delr) ;
            end
            
            if (j==1)
                term2= (V(i,j+1)-V(i,t_n))/(2*delth);
            elseif (j==t_n)
                term2= (V(i,1)-V(i,j-1))/(2*delth);
            else
                term2= (V(i,j+1)-V(i,j-1))/(2*delth);
            end
                term3 = U(i,j);
            deldotU(i,j) =  term1 + (1/r(i))*(term2+term3);
        end
    end
  
    
    b=zeros(size(XX));
    
    for i=2:r_n-1 %1 to r_n
        for j=1:t_n
            b(i,j)=(1/dt)*deldotU(i,j);
        end
    end
    
%     
%     
%  %%%   b(r_n-3,t_n-3)=0;
%     
%     
%     % Biconjugate Gradient Stabilized method %%%%%%%%%%%%%%%%%%%
%     %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     % Vectors for BICGSTAB
%     v0=zeros(size(XX));
%     d0=zeros(size(XX));
%     dd0=zeros(size(XX));
%     x0=zeros(size(XX));%zeros(size(XX));
%     %x0(10,10)=1;
%     % Scalars
%     p0=v0;
%     alpha =1;
%     w0=1;
%     rho0=1;
%     
%     iter = 0;
%     
%     Ax0=ComputeAx_collocated_pn(x0,rr,tt, delr, delth, r_n, t_n,t,Pneumann);
%     
%     for i=2:r_n-1
%         for j=1:t_n
%             dd0(i,j)=b(i,j)-Ax0(i,j);
%             d0(i,j)=dd0(i,j);
%         end
%     end
%     %dd0=b-Ax0;
%     
%     %d0=dd0;
%     
%     
%     Tol=1e-15;
%     Res=100;
%     
%     while (Res>Tol && iter<5000),
%         rho_i=innerp(dd0,d0,r_n,t_n);
%         if (rho_i==0) break, end
%         if (iter>0)
%             beta = (rho_i/rho0)*(alpha/w0);
%             p_i = d0 + beta.*(p0 - w0.*v0);
%             
%         else
%             p_i=d0;
%             
%         end
%         v_i = ComputeAx_collocated_pn(p_i,rr,tt,delr,delth,r_n,t_n,t,Pneumann);
%         
%         alpha = rho_i/innerp(dd0,v_i,r_n,t_n);
%         s = d0 - alpha.*v_i;
%         
%         t=ComputeAx_collocated_pn(s,rr,tt,delr,delth,r_n,t_n,t,Pneumann);
%         
%         w_i=innerp(t,s,r_n,t_n)/innerp(t,t,r_n,t_n);
%         x_i=x0+alpha.*p_i + w_i.*s;
%         
%         d_i=s - w_i.*t;
%         
%         P= x_i;
%         %%%%%%%%
%         d0=d_i;
%         Res = innerp((x_i - x0 ), (x_i -x0),r_n,t_n);
%         Res = sqrt ( Res / ((r_n-2)*(t_n)));
%         x0=x_i;
%         p0=p_i;
%         w0=w_i;
%         v0=v_i;
%         rho0=rho_i;
%         
%         iter=iter+1;
%         
%     
%     end
%     fprintf('\n Pressure converged at Iter %d, Residual is %0.3g at time step %d', iter, Res, k);
    % ------------------------------------------------------------
    
    %b=ones(t_n*r_n,1);
    
    
    P = BICGSTABIter_p(P,r,r_n,t_n,delr,delth,b);

    % Making the plot look good by accounting for Neumann boundary
    for j=1:t_n

        P(1,j)=P(2,j);
       P(r_n,j) =P(r_n-1,j);
    end
    
    
    for i=2:r_n-1
        for j=1:t_n
            
                G1(i,j) = (P(i+1,j)-P(i-1,j))/(2*delr);
            
            if (j==1)
                G2(i,j) =(1/r(i))*(P(i,j+1)-P(i,t_n))/(2*delth);
            elseif (j==t_n)
                G2(i,j) =(1/r(i))*(P(i,1)-P(i,j-1))/(2*delth);
            else
                G2(i,j) =(1/r(i))*(P(i,j+1)-P(i,j-1))/(2*delth);
            end
 
            U(i,j)= U(i,j) - dt*G1(i,j);
            V(i,j)= V(i,j) - dt*G2(i,j);
        end
    end
   
       
    for i=1:r_n
        for j=1:t_n+1
            if j==t_n+1
                Uplot(i,j) = U(i,1);
                Vplot(i,j) = V(i,1);
                Pplot(i,j) = P(i,1);
            else
                Uplot(i,j)=U(i,j); Vplot(i,j) = V(i,j);
                Pplot(i,j) = P(i,j);
            end
            
        end
    end
    
     
    sum=0;
    for j=1:t_n
        
        sum= sum+U(r_n,j);
    end
    fprintf('\n Sum after Pressure correction = %f \n',sum);
    
    if (mod(k,writenum)==0)
        %Velocities in XY coordinates for the contours.
        for i=1:r_n
            for j=1:t_n+1
                Umag(i,j)=sqrt(Uplot(i,j)^2 + Vplot(i,j)^2);
                if j==t_n+1
                Ux(i,j)=Uplot(i,j)*cos(tt(i,1)) - Vplot(i,j)*sin(tt(i,1));
                Uy(i,j)=Uplot(i,j)*sin(tt(i,1))+ Vplot(i,j)*cos(tt(i,1));
                else
                    Ux(i,j)=Uplot(i,j)*cos(tt(i,j)) - Vplot(i,j)*sin(tt(i,j));
                Uy(i,j)=Uplot(i,j)*sin(tt(i,j))+ Vplot(i,j)*cos(tt(i,j));
                end
            end
        end
        
        % writing the solutions
        if k==0
            fid2 = fopen('sol_init.xyz','w');
        elseif k==1000
            fid2 = fopen('sol_1000.xyz','w');
        elseif k==2500
            fid2 = fopen('sol_2500.xyz','w');
        elseif k==5000
            fid2 = fopen('sol_5000.xyz','w');
        elseif k==8000
            fid2 = fopen('sol_8000.xyz','w');
        elseif k==11000
            fid2 = fopen('sol_11000.xyz','w');
        elseif k==14000
            fid2 = fopen('sol_14000.xyz','w');
        elseif k==17000
            fid2 = fopen('sol_17000.xyz','w');
        elseif k==20000
            fid2 = fopen('sol_20000.xyz','w');
        elseif k==25000
            fid2 = fopen('sol_25000.xyz','w');
        elseif k==30000
            fid2 = fopen('sol_30000.xyz','w');
        elseif k==35000
            fid2 = fopen('sol_35000.xyz','w');
        elseif k==40000
            fid2 = fopen('sol_40000.xyz','w');
        else
            fid2 = fopen('sol.xyz','w');
        end
        fprintf(fid2, '%d \t %d \t %d\n', r_n, t_n+1, 6);
        % Pressure
        fprintf(fid2,'%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n',Pplot);
        fprintf(fid2,'%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n',Uplot);
        fprintf(fid2,'%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n',Vplot);
        fprintf(fid2,'%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n',Ux);
        fprintf(fid2,'%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n',Uy);
        fprintf(fid2,'%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n',Umag);
        % U
        % V
        fclose(fid2);
        fprintf('\n Wrote solution to file \n');
    end
    
    
end %time###########################################################3

end
