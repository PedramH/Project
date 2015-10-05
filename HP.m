%Version 1.0   (High Performance)

clear all
clc

%Global variables 
dt = 0.001;
etta = 9;
detta = 0.05;
%----------------
C = 0.5;
phi = 0.2;		        %Solid volume fraction of the nanofluid 

%------------------Base Fluid-----------------------

pf = 997.1 ;			%Density 
Kf  = 0.613 ;			%Thermal conductivity 
pCpf = pf * 4179;		%Heat capacitance 
%Vf = 0.8926e-06;		%Kinematic viscosity (m2/s)
Pr = 1 ;                %Prandtl number 5.78              
%-----------------NanoParticle----------------------

NanoParticle = 'Cu' ;		%Valid names : Cu  Al2O3  TiO2 

switch (NanoParticle)

   case 'Cu'
       ps =  8933 ;		    %Density  
       Ks = 400;			%Thermal conductivity 
       pCps = ps * 385;		%Heat capacitance
   case 'Al2O3'
       ps =  3970 ;		      
       Ks = 40;			
       pCps = ps * 765;		
   case 'TiO2'
       ps =  4250 ;		     
       Ks = 8.9538;			 
       pCps = ps * 686.2;		
   otherwise
	   error('Wrong NanoParticle name.')

end

%------------------NanoFluid------------------------

Knf =( ((Ks+2*Kf)-2*phi*(Kf-Ks))/((Ks+2*Kf)+phi*(Kf-Ks)) ).*Kf;					%Thermal conductivity 
pCpnf = (1-phi).*(pCpf) + phi .*(pCps) ; 										%Heat capacitance
anf = Knf / pCpnf ; 															%Effective thermal diffusivity 

%---------------------------------------------------

e1 = 1/( ((1-phi).^(2.5)).*((1-phi) + phi.*(ps/pf)) );
e2 = (Knf/Kf)/ ( (1-phi)+phi.*(pCps/pCpf) );   
%Pr = Vf/anf; 




tic
for kdt = 1:1000   %main loop  -- 
   
    
    n = etta/detta;
    % initial F, W, G values
    F = zeros(n,2);
	G = zeros(n,2);
    W = zeros(n,2);
	
    f = zeros(n,1);
	g = zeros(n,1);
    
    nt = 1000000000000;  %No. of iterations
   
    % K initialization
    k1F = zeros(n,1); k2F = zeros(n,1); k3F = zeros(n,1); k4F = zeros(n,1);
	k1G = zeros(n,1); k2G = zeros(n,1); k3G = zeros(n,1); k4G = zeros(n,1);
    k1W = zeros(n,1); k2W = zeros(n,1); k3W = zeros(n,1); k4W = zeros(n,1);
    
    % B.C   
    F(1,:) = 1; G(1,:) = C; W(1,:) = 1;
    F(n,:) = 0; G(n,:) = 0; W(n,:) = 0;
    f(1) = 0; g(1) = 0;
    
   
   for k=1:nt    %Time Loop
       
	  
       
           i= 2:n-1;  
		   
		   %K1
           
		   f(i)=f(i-1)+F(i,1).*detta; 
		   g(i)=g(i-1)+G(i,1).*detta;
		   
		   
           %-------------------------F----------------------------------
           k1F(i) = dt.* (    e1.*( (F(i+1,1) - 2.*F(i,1) + F(i-1,1) )./(detta.^2) )+...
							(f(i)+g(i)).*( (F(i,1)-F(i-1,1))./detta ) -...
							F(i,1).^2         );	
           
                        
           %-------------------------G----------------------------------
           k1G(i) = dt.* (    e1.*( (G(i+1,1) - 2.*G(i,1) + G(i-1,1) )./(detta.^2) )+...
							(f(i)+g(i)).*( (G(i,1)-G(i-1,1))./detta ) -...
							G(i,1).^2         )	;
           
		   %-------------------------W----------------------------------
           k1W(i) = dt.* (   (e2./Pr).* ( (W(i+1,1) - 2.*W(i,1) + W(i-1,1) )./(detta.^2) ) +...
							(f(i)+g(i)).*( (W(i,1)- W(i-1,1))./detta ));
						
		  
           
       
     
       
      
       
       
       
       
           %K2
           
           f(i)=f(i-1)+(F(i,1)+k1F(i).*0.5).*detta;
           g(i)=g(i-1)+(G(i,1)+k1G(i).*0.5).*detta;
		  
		   %-------------------------F----------------------------------
           k2F(i) = dt.* (    e1.*( (F(i+1,1)+k1F(i+1).*0.5 - 2.*(F(i,1)+k1F(i).*0.5) + F(i-1,1)+k1F(i-1).*0.5 )./(detta.^2) )+...
							(f(i)+g(i)).*( (F(i,1)+k1F(i).*0.5 - (F(i-1,1)+k1F(i-1).*0.5))./detta ) -...
							(F(i,1)+k1F(i).*0.5) .^ 2         )	;
           
           %-------------------------G----------------------------------
           k2G(i) = dt.* (    e1.*( (G(i+1,1)+k1G(i+1).*0.5 - 2.*(G(i,1)+k1G(i).*0.5) + G(i-1,1)+k1G(i-1).*0.5 )./(detta.^2) )+...
							(f(i)+g(i)).*( (G(i,1)+k1G(i).*0.5 - (G(i-1,1)+k1G(i-1).*0.5))./detta ) -...
							(G(i,1)+k1G(i).*0.5) .^ 2         )	;
		   
		   %-------------------------W----------------------------------
           k2W(i) = dt.* (   (e2./Pr).* ( (W(i+1,1)+k1W(i+1).*0.5 - 2.*(W(i,1)+k1W(i).*0.5) + W(i-1,1)+k1W(i-1).*0.5 )./(detta.^2) ) +...
							(f(i)+g(i)).*( (W(i,1)+k1W(i).*0.5 - (W(i-1,1)+k1W(i-1).*0.5))./detta )  );
      
       
       
       
       
       
       
       
       
    
           %K3
           
           f(i)=f(i-1)+(F(i,1)+k2F(i).*0.5).*detta;
           g(i)=g(i-1)+(G(i,1)+k2G(i).*0.5).*detta;
		  
		   %-------------------------F----------------------------------
           k3F(i) = dt.* (    e1.*( (F(i+1,1)+k2F(i+1).*0.5 - 2.*(F(i,1)+k2F(i).*0.5) + F(i-1,1)+k2F(i-1).*0.5 )./(detta.^2) )+...
							(f(i)+g(i)).*( (F(i,1)+k2F(i).*0.5 - (F(i-1,1)+k2F(i-1).*0.5))./detta ) -...
							(F(i,1)+k2F(i).*0.5) .^ 2         )	;
           
           %-------------------------G----------------------------------
           k3G(i) = dt.* (    e1.*( (G(i+1,1)+k2G(i+1).*0.5 - 2.*(G(i,1)+k2G(i).*0.5) + G(i-1,1)+k2G(i-1).*0.5 )./(detta.^2) )+...
							(f(i)+g(i)).*( (G(i,1)+k2G(i).*0.5 - (G(i-1,1)+k2G(i-1).*0.5))./detta ) -...
							(G(i,1)+k2G(i).*0.5) .^ 2         )	;
		   
		   %-------------------------W----------------------------------
           k3W(i) = dt.* (   (e2./Pr).* ( (W(i+1,1)+k2W(i+1).*0.5 - 2.*(W(i,1)+k2W(i).*0.5) + W(i-1,1)+k2W(i-1).*0.5 )./(detta.^2) ) +...
							(f(i)+g(i)).*( (W(i,1)+k2W(i).*0.5 - (W(i-1,1)+k2W(i-1).*0.5))./detta )  );
       
       
       
       
       
       
       

       
           %K4
           
           f(i)=f(i-1)+(F(i,1)+k3F(i).*1).*detta;
           g(i)=g(i-1)+(G(i,1)+k3G(i).*1).*detta;
		  
		   %-------------------------F----------------------------------
           k4F(i) = dt.* (    e1.*( (F(i+1,1)+k3F(i+1).*1 - 2.*(F(i,1)+k3F(i).*1) + F(i-1,1)+k3F(i-1).*1 )./(detta.^2) )+...
							(f(i)+g(i)).*( (F(i,1)+k3F(i).*1 - (F(i-1,1)+k3F(i-1).*1))./detta ) -...
							(F(i,1)+k3F(i).*1) .^ 2         )	;
           
           %-------------------------G----------------------------------
           k4G(i) = dt.* (    e1.*( (G(i+1,1)+k3G(i+1).*1 - 2.*(G(i,1)+k3G(i).*1) + G(i-1,1)+k3G(i-1).*1 )./(detta.^2) )+...
							(f(i)+g(i)).*( (G(i,1)+k3G(i).*1 - (G(i-1,1)+k3G(i-1).*1))./detta ) -...
							(G(i,1)+k3G(i).*1) .^ 2         )	;
		   
		   %-------------------------W----------------------------------
           k4W(i) = dt.* (   (e2./Pr).* ( (W(i+1,1)+k3W(i+1).*1 - 2.*(W(i,1)+k3W(i).*1) + W(i-1,1)+k3W(i-1).*1 )./(detta.^2) ) +...
							(f(i)+g(i)).*( (W(i,1)+k3W(i).*1 - (W(i-1,1)+k3W(i-1).*1))./detta )  );
       
       
       
       
       
	       %Caclulate F,G,W in t=t+dt
       
           F(i,2)=F(i,1)+(1./6).*(k1F(i)+2.*k2F(i)+2.*k3F(i)+k4F(i));
		   G(i,2)=G(i,1)+(1./6).*(k1G(i)+2.*k2G(i)+2.*k3G(i)+k4G(i));
           W(i,2)=W(i,1)+(1./6).*(k1W(i)+2.*k2W(i)+2.*k3W(i)+k4W(i));
       
       
       % B.C   
       % F(1,2) = 1; G(1,2) = C; W(1,2) = 1;
       % F(n,2) = 0; G(n,2) = 0; W(n,2) = 0;
       
       NomF = norm((F(:,1)-F(:,2)))/norm(F(:,1))
       if  NomF < 1e-7 && k > 2000
           disp('Checkpoint 1')
           break
       end
       
       
       %Check for divergence
       if isnan(NomF)
           dt = dt - 0.2*dt
           clear F G W f g
           break
       end
	   
	   
	   F(:,1)=F(:,2);
	   G(:,1)=G(:,2);
	   W(:,1)=W(:,2);
	   
           
       
   end
    
   if  NomF < 1e-7 && k > 2000
       disp('Checkpoint 2')
       kdt
       break
   end
   
end
time = toc
x=0:detta:etta-detta;
%----------------Save Output Data---------------------------
No = strrep(num2str(phi), '.', '_');
name = sprintf('Data\\%s\\phi%s.mat',NanoParticle,No);
save(name)
