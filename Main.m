%version 0.1

clear all
clc

%Global variables 
dt = 0.0001;
etta = 5;
detta = 0.05;
%----------------
C = 0.5;
e1 = 1;
e2 = 1;
Pr = 1;




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
    F(1,1) = 1; G(1,1) = C; W(1,1) = 1;
    F(n,1) = 0; G(n,1) = 0; W(n,1) = 0;
    f(1) = 0; g(1) = 0;
    
   
   for k=1:nt    %Time Loop
       
	  
       
       for i= 2:n-1  %K1
           
           f(i)=f(i-1)+F(i,1)*detta; 
		   g(i)=g(i-1)+G(i,1)*detta;
		   
		   
           %-------------------------F----------------------------------
           k1F(i) = dt* (    e1*( (F(i+1,1) - 2*F(i,1) + F(i-1,1) )/(detta.^2) )+...
							(f(i)+g(i))*( (F(i,1)-F(i-1,1))/detta ) -...
							F(i,1).^2         );	
           
                        
           %-------------------------G----------------------------------
           k1G(i) = dt* (    e1*( (G(i+1,1) - 2*G(i,1) + G(i-1,1) )/(detta.^2) )+...
							(f(i)+g(i))*( (G(i,1)-G(i-1,1))/detta ) -...
							G(i,1).^2         )	;
           
		   %-------------------------W----------------------------------
           k1W(i) = dt* (   (e2/Pr)* ( (W(i+1,1) - 2*W(i,1) + W(i-1,1) )/(detta.^2) ) +...
							(f(i)+g(i))*( (W(i,1)- W(i-1,1))/detta ));
						
		  
           
       end
     
       
      
       
       
       
       
       for i= 2:n-1  %K2
           % TODO
           f(i)=f(i-1)+(F(i,1)+k1F(i)*0.5)*detta;
           g(i)=g(i-1)+(G(i,1)+k1G(i)*0.5)*detta;
		  
		   %-------------------------F----------------------------------
           k2F(i) = dt* (    e1*( (F(i+1,1)+k1F(i+1)*0.5 - 2*(F(i,1)+k1F(i)*0.5) + F(i-1,1)+k1F(i-1)*0.5 )/(detta.^2) )+...
							(f(i)+g(i))*( (F(i,1)+k1F(i)*0.5 - (F(i-1,1)+k1F(i-1)*0.5))/detta ) -...
							(F(i,1)+k1F(i)*0.5) .^ 2         )	;
           
           %-------------------------G----------------------------------
           k2G(i) = dt* (    e1*( (G(i+1,1)+k1G(i+1)*0.5 - 2*(G(i,1)+k1G(i)*0.5) + G(i-1,1)+k1G(i-1)*0.5 )/(detta.^2) )+...
							(f(i)+g(i))*( (G(i,1)+k1G(i)*0.5 - (G(i-1,1)+k1G(i-1)*0.5))/detta ) -...
							(G(i,1)+k1G(i)*0.5) .^ 2         )	;
		   
		   %-------------------------W----------------------------------
           k2W(i) = dt* (   (e2/Pr)* ( (W(i+1,1)+k1W(i+1)*0.5 - 2*(W(i,1)+k1W(i)*0.5) + W(i-1,1)+k1W(i-1)*0.5 )/(detta.^2) ) +...
							(f(i)+g(i))*( (W(i,1)+k1W(i)*0.5 - (W(i-1,1)+k1W(i-1)*0.5))/detta )  );
       end
       
       
       
       
       
       
       
       
    
       for i= 2:n-1  %K3
           % TODO
           f(i)=f(i-1)+(F(i,1)+k2F(i)*0.5)*detta;
           g(i)=g(i-1)+(G(i,1)+k2G(i)*0.5)*detta;
		  
		   %-------------------------F----------------------------------
           k3F(i) = dt* (    e1*( (F(i+1,1)+k2F(i+1)*0.5 - 2*(F(i,1)+k2F(i)*0.5) + F(i-1,1)+k2F(i-1)*0.5 )/(detta.^2) )+...
							(f(i)+g(i))*( (F(i,1)+k2F(i)*0.5 - (F(i-1,1)+k2F(i-1)*0.5))/detta ) -...
							(F(i,1)+k2F(i)*0.5) .^ 2         )	;
           
           %-------------------------G----------------------------------
           k3G(i) = dt* (    e1*( (G(i+1,1)+k2G(i+1)*0.5 - 2*(G(i,1)+k2G(i)*0.5) + G(i-1,1)+k2G(i-1)*0.5 )/(detta.^2) )+...
							(f(i)+g(i))*( (G(i,1)+k2G(i)*0.5 - (G(i-1,1)+k2G(i-1)*0.5))/detta ) -...
							(G(i,1)+k2G(i)*0.5) .^ 2         )	;
		   
		   %-------------------------W----------------------------------
           k3W(i) = dt* (   (e2/Pr)* ( (W(i+1,1)+k2W(i+1)*0.5 - 2*(W(i,1)+k2W(i)*0.5) + W(i-1,1)+k2W(i-1)*0.5 )/(detta.^2) ) +...
							(f(i)+g(i))*( (W(i,1)+k2W(i)*0.5 - (W(i-1,1)+k2W(i-1)*0.5))/detta )  );
       end
       
       
       
       
       
       

       
       for i= 2:n-1  %K4
           % TODO
           f(i)=f(i-1)+(F(i,1)+k3F(i)*1)*detta;
           g(i)=g(i-1)+(G(i,1)+k3G(i)*1)*detta;
		  
		   %-------------------------F----------------------------------
           k4F(i) = dt* (    e1*( (F(i+1,1)+k3F(i+1)*1 - 2*(F(i,1)+k3F(i)*1) + F(i-1,1)+k3F(i-1)*1 )/(detta.^2) )+...
							(f(i)+g(i))*( (F(i,1)+k3F(i)*1 - (F(i-1,1)+k3F(i-1)*1))/detta ) -...
							(F(i,1)+k3F(i)*1) .^ 2         )	;
           
           %-------------------------G----------------------------------
           k4G(i) = dt* (    e1*( (G(i+1,1)+k3G(i+1)*1 - 2*(G(i,1)+k3G(i)*1) + G(i-1,1)+k3G(i-1)*1 )/(detta.^2) )+...
							(f(i)+g(i))*( (G(i,1)+k3G(i)*1 - (G(i-1,1)+k3G(i-1)*1))/detta ) -...
							(G(i,1)+k3G(i)*1) .^ 2         )	;
		   
		   %-------------------------W----------------------------------
           k4W(i) = dt* (   (e2/Pr)* ( (W(i+1,1)+k3W(i+1)*1 - 2*(W(i,1)+k3W(i)*1) + W(i-1,1)+k3W(i-1)*1 )/(detta.^2) ) +...
							(f(i)+g(i))*( (W(i,1)+k3W(i)*1 - (W(i-1,1)+k3W(i-1)*1))/detta )  );
       end
       
       
       
       
	   %Caclulate F,W in t=t+dt
       for i=2:n-1
           F(i,2)=F(i,1)+(1/6)*(k1F(i)+2*k2F(i)+2*k3F(i)+k4F(i));
		   G(i,2)=G(i,1)+(1/6)*(k1G(i)+2*k2G(i)+2*k3G(i)+k4G(i));
           W(i,2)=W(i,1)+(1/6)*(k1W(i)+2*k2W(i)+2*k3W(i)+k4W(i));
       end
       
       % B.C   
       F(1,2) = 1; G(1,2) = C; W(1,2) = 1;
       F(n,2) = 0; G(n,2) = 0; W(n,2) = 0;
       
       NomF = norm((F(:,1)-F(:,2)))/norm(F(:,1))
       if  k > 2000 && NomF < 1e-7
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
    
   if  k > 2000 && NomF < 1e-7
       disp('Checkpoint 2')
       kdt
       break
   end
   
end
time = toc
save('Test1.mat')
x=0:detta:etta-detta;
% plot(x,F(:,k))
% plot(x,W(:,k))