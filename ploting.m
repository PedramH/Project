%Plot

clear 
clc 


name = sprintf('Data\\AllData.mat');
load(name)
%%% AllData.mat
% Fs 3*n*FileNo     [Nano Particle , eta , phi]  <-- Steady state
% Fu 3*n*j*FileNo   [Nano Particle , eta , time , phi] <--  unsteady 
% 

in = 1 ;  NanoParticle = NanoParticles{in};  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Plot F,G & theta with diffrent phi(steady state)      %%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(0)   %RUN? 
titlename = sprintf(' NanoParticle : %s' ,NanoParticle);
h(1)=figure(1); %F
plot(x,Fs(in,:,1),...
    x,Fs(in,:,6),...
    x,Fs(in,:,11),...
    x,Fs(in,:,21),...
    'LineWidth',2);
title(titlename,'FontSize',20)
xlabel('\eta','FontSize',20)
ylabel('F','FontSize',20);
h_legend(1) =legend('\phi=0','\phi=0.05','\phi=0.1','\phi=0.2');
set(h_legend(1),'FontSize',20);
set(gca,'fontsize',15)
name = sprintf('Data\\%s\\figs\\F',NanoParticle);
print(name,'-djpeg','-r500')
%-----------------------------------------------------------
h(2)=figure(2); %G
plot(x,Gs(in,:,1),...
    x,Gs(in,:,6),...
    x,Gs(in,:,11),...
    x,Gs(in,:,21),...
    'LineWidth',2);
title(titlename,'FontSize',20)
xlabel('\eta','FontSize',20)
ylabel('G','FontSize',20);
h_legend(2)=legend('\phi=0','\phi=0.05','\phi=0.1','\phi=0.2');
set(h_legend(2),'FontSize',20);
name = sprintf('Data\\%s\\figs\\G',NanoParticle);
print(name,'-djpeg','-r500')
%------------------------------------------------------------
h(3)=figure(3); %W
plot(x,Ws(in,:,1),...
    x,Ws(in,:,6),...
    x,Ws(in,:,11),...
    x,Ws(in,:,21),...
    'LineWidth',2);
title(titlename,'FontSize',20)
xlabel('\eta','FontSize',20)
ylabel('\theta','FontSize',20);
h_legend(3)=legend('\phi=0','\phi=0.05','\phi=0.1',...
    '\phi=0.2');

set(h_legend(3),'FontSize',20);
name = sprintf('Data\\%s\\figs\\theta',NanoParticle);
print(name,'-djpeg','-r500')

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Plot F,G & theta with diffrent time and phi           %%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (1) % RUN? 
 clear
 name = sprintf('Data\\AllData_Old.mat');
 load(name)
 phi = 0;
 in = 3 ;  NanoParticle = NanoParticles{in};  
 for i=1:FileNo
    
    
    %----------------------F---------------------%
    figure(1);
    phi_str = num2str(phi);
    titlename = sprintf(' NanoParticle : %s \n\\phi = %s',NanoParticle,phi_str);
    
    plot(x,Fu(in,:,1,i),...
        x,Fu(in,:,3,i),...
        x,Fu(in,:,8,i),...
        x,Fu(in,:,12,i),...
        x,Fu(in,:,15,i),...
        x,Fs(in,:,i),'LineWidth',2)
    title(titlename,'FontSize',20);
    xlabel('\eta','FontSize',20);
    ylabel('F','FontSize',20);
    %legend('t*= 1','2','3','4','5','Steady')
    
    % File name 
    No = strrep(num2str(phi), '.', '_');
    name = sprintf('Data\\%s\\figs\\F\\phi%s',NanoParticle,No);
    
    % Create arrow
    annotation('arrow',[0.16 0.39],...
    [0.18 0.32]);
   
    % Create Text 
    text('String','t* = 0.02, 0.06, 0.16, 0.6, 1.2, Steady',...
    'Position',[1.72235023041475 0.258982035928144 17.3205080756888],...
    'FontSize',14);
    % Save figure
    print(name,'-djpeg','-r300')
    
    %----------------------G---------------------%
    figure(2);
    phi_str = num2str(phi);
    titlename = sprintf(' NanoParticle : %s \n\\phi = %s',NanoParticle,phi_str);
    
    plot(x,Gu(in,:,1,i),...
        x,Gu(in,:,3,i),...
        x,Gu(in,:,8,i),...
        x,Gu(in,:,12,i),...
        x,Gu(in,:,15,i),...
        x,Gs(in,:,i),'LineWidth',2)
    title(titlename,'FontSize',20);
    xlabel('\eta','FontSize',20);
    ylabel('G','FontSize',20);
    %legend('t*= 1','2','3','4','5','Steady')
    
    % File name 
    name = sprintf('Data\\%s\\figs\\G\\phi%s',NanoParticle,No);
    
    % Create arrow
    annotation('arrow',[0.16 0.39],...
    [0.18 0.32]);
   
    % Create Text 
    text('String','t* = 0.02, 0.06, 0.16, 0.6, 1.2, Steady',...
    'Position',[1.73387096774194 0.135479041916168 17.3205080756888],...
    'FontSize',14);
    % Save figure
    print(name,'-djpeg','-r300')
    
    
    %----------------------theta---------------------%
    figure(3);
    phi_str = num2str(phi);
    titlename = sprintf(' NanoParticle : %s \n\\phi = %s',NanoParticle,phi_str);
    
    plot(x,Wu(in,:,1,i),...
        x,Wu(in,:,3,i),...
        x,Wu(in,:,8,i),...
        x,Wu(in,:,12,i),...
        x,Wu(in,:,15,i),...
        x,Ws(in,:,i),'LineWidth',2)
    title(titlename,'FontSize',20);
    xlabel('\eta','FontSize',20);
    ylabel('\theta','FontSize',20);
    %legend('t*= 1','2','3','4','5','Steady')
    
    % File name 
    name = sprintf('Data\\%s\\figs\\theta\\phi%s',NanoParticle,No);
    
    % Create arrow
    annotation('arrow',[0.14 0.30],...
    [0.17 0.30]);
   
    % Create Text 
    text('String','t* = 0.02, 0.06, 0.16, 0.6, 1.2, Steady',...
    'Position',[1.16596638655462 0.238911290322581 17.3205080756888],...
    'FontSize',14);
    % Save figure
    print(name,'-djpeg','-r300')   
    
    
    phi = phi + dphi;
 end
 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Plot Nu and Cf with respect to phi                   %%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if (0)   %RUN? 
    
    Cfx=zeros(3,FileNo); Cfy = zeros(3,FileNo);
    xphi=0:dphi:MAXphi-dphi;
    
    for in=1:3
        phi = 0;
        for i=1:FileNo
            
            dF = (Fs(in,2,i)-Fs(in,1,i))/detta;
            dG = (Gs(in,2,i)-Gs(in,1,i))/detta;
            dW = (Ws(in,2,i)-Ws(in,1,i))/detta;
            
            
            Cfx(in,i) = -dF/(1-phi).^(2.5);
            Cfy(in,i) = -dG/( ((1-phi).^(2.5))*C.^(3/2) );
            
            Nux(in,i) = -dW * (Knf(in,i)/Kf);
            Nuy(in,i) = -dW * (Knf(in,i)/(Kf*C.^(1/2)));
        
            phi = phi + dphi; 
        end
    end
    figure(1)
    hold on
    
    xlabel('\phi','FontSize',20)
    ylabel('C_{f}','FontSize',20)
    
    h1 = plot(xphi,Cfx(1,:),'k',xphi,Cfy(1,:),'--k','LineWidth',2);
    
    
    h2 = plot(xphi,Cfx(2,:),'b',xphi,Cfy(2,:),'--b','LineWidth',2);
    

    h3 = plot(xphi,Cfx(3,:),'r',xphi,Cfy(3,:),'--r','LineWidth',2);
    
    text('String','Cu',...
    'Position',[0.0525604838709678 2.78434297868618 17.3205080756888],...
    'FontSize',18,'Color','k');
    text('String','Al2O3',...
    'Position',[0.0712298387096774 2.78434297868618 17.3205080756888],...
    'FontSize',18,'Color','b');
    text('String','TiO2',...
    'Position',[0.096754032258065 2.78434297868618 17.3205080756888],...
    'FontSize',18,'Color','r');


    h1_legend=legend(h1,'Re_{x}^{1/2}C_{fx}','Re_{y}^{1/2}C_{fy}','Location','best');
    set(h1_legend,'FontSize',15)
    hold off
    name = sprintf('Data\\comp\\Cf');
    print(name,'-djpeg','-r300')
    
    figure(2)
    plot(xphi,Nux(1,:),'k',xphi,Nuy(1,:),'--k','LineWidth',2);
    xlabel('\phi','FontSize',20)
    ylabel('Nu','FontSize',20)
    name = sprintf('Data\\comp\\Nu');
    print(name,'-djpeg','-r300')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Plot Nu and Cf with respect to time                   %%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if (0)   % RUN? 
    j=300;
    t=10:1:j-1;
    t = t*20*dt;
    Cfx=zeros(3,290); Cfy = zeros(3,290);Nux=zeros(3,290);Nuy=zeros(3,290);
    
%     in = 1 ; % Cu 
%     NanoParticle = NanoParticles{in};
    
    phi = 0;
    for i=1:FileNo  %Change this 
        for in = 1:3
            for k = 1:j-10    
                dF = (Fu(in,2,k+10,i)-Fu(in,1,k+10,i))/detta;
                dG = (Gu(in,2,k+10,i)-Gu(in,1,k+10,i))/detta;
                dW = (Wu(in,2,k+10,i)-Wu(in,1,k+10,i))/detta;
                Cfx(in,k) = -dF/(1-phi).^(2.5);
                Cfy(in,k) = -dG/( ((1-phi).^(2.5))*C.^(3/2) );
            
                Nux(in,k) = -dW * (Knf(in,i)/Kf);
                Nuy(in,k) = -dW * (Knf(in,i)/(Kf*C.^(1/2)));
            end
        end
        
        
        
        
        
        
        % ------------Cf---------------
        figure(1)
        plot(t,Cfx(1,:),'r',t,Cfx(2,:),'g',...
            t,Cfx(3,:),'b',t,Cfy(1,:),'--r',...
            t,Cfy(2,:),'--g',t,Cfy(3,:),'--b',...
            'LineWidth',2)
        
        legend({'Cu','Al_{2}O_{3}','TiO_{2}'}...
            ,'FontSize',15);
        xlabel('t*','FontSize',20)
        ylabel('C_{f}','FontSize',20)
%         phi_str = num2str(phi);
%         titlename = sprintf(' NanoParticle : %s \n\\phi = %s',NanoParticle,phi_str);
%         title('Re_{x}^{1/2}C_{f}','FontSize',20)
        xlim([0.2,inf])
        
        % save figure 
        No = strrep(num2str(phi), '.', '_');
        name = sprintf('Data\\All\\Cf\\phi%s',No);
        print(name,'-djpeg','-r300')
        
        % ------------Nu---------------
        figure(2)
        plot(t,Nux(1,:),'r',t,Nux(2,:),'g',...
            t,Nux(3,:),'b',t,Nuy(1,:),'--r',...
            t,Nuy(2,:),'--g',t,Nuy(3,:),'--b',...
            'LineWidth',2)
        
        legend({'Cu','Al_{2}O_{3}','TiO_{2}'}...
            ,'FontSize',15);
        xlabel('t*','FontSize',20)
        ylabel('Nu','FontSize',20)
%         phi_str = num2str(phi);
%         titlename = sprintf(' NanoParticle : %s \n\\phi = %s',NanoParticle,phi_str);
%         title(titlename,'FontSize',20)
        %ylim([-0.2,5])
        xlim([0.2,inf])
        % save figure 
        No = strrep(num2str(phi), '.', '_');
        name = sprintf('Data\\All\\Nu\\phi%s',No);
        print(name,'-djpeg','-r300')
        
        
        phi = phi + dphi; 
    end
    
 end
    