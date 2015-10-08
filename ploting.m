%Plot

clear 
clc 

NanoParticle = 'Cu';
name = sprintf('Data\\%s\\AllData.mat',NanoParticle);
load(name)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Plot F,G & theta with diffrent phi(steady state)      %%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(0)   %RUN? 

h(1)=figure(1); %F
plot(x,Fs(:,1),x,Fs(:,7),x,Fs(:,14)...
    ,x,Fs(:,21));
xlabel('\eta')
ylabel('F');
legend('\phi=0','\phi=0.06','\phi=0.13','\phi=0.2')

h(2)=figure(2); %G
plot(x,Gs(:,1),x,Gs(:,7),x,Gs(:,14)...
    ,x,Gs(:,21));
xlabel('\eta')
ylabel('G');
legend('\phi=0','\phi=0.06','\phi=0.13','\phi=0.2')

h(3)=figure(3); %W
plot(x,Ws(:,1),x,Ws(:,3),x,Ws(:,7),x,Ws(:,10),x,Ws(:,14)...
    ,x,Ws(:,21));
xlabel('\eta')
ylabel('\theta');
legend('\phi=0','\phi=0.02','\phi=0.06','\phi=0.09','\phi=0.13','\phi=0.2')

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Plot F,G & theta with diffrent time and phi           %%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (1) % RUN? 

 phi = 0;
 for i=1:FileNo
    
    %----------------------F---------------------%
    figure(1);
    phi_str = num2str(phi);
    titlename = sprintf(' NanoParticle : %s \n\\phi = %s',NanoParticle,phi_str);
    
    plot(x,Fu(:,1,i),...
        x,Fu(:,3,i),...
        x,Fu(:,8,i),...
        x,Fu(:,12,i),...
        x,Fu(:,15,i),...
        x,Fs(:,i),'LineWidth',1.2)
    title(titlename);
    xlabel('\eta');
    ylabel('F');
    %legend('t*= 1','2','3','4','5','Steady')
    
    % File name 
    No = strrep(num2str(phi), '.', '_');
    name = sprintf('Data\\%s\\figs\\F\\phi%s',NanoParticle,No);
    
    % Create arrow
    annotation('arrow',[0.16 0.39],...
    [0.18 0.32]);
   
    % Create Text 
    text('String','t*',...
    'Position',[1.72235023041475 0.258982035928144 17.3205080756888],...
    'FontSize',18);
    % Save figure
    print(name,'-djpeg','-r300')
    
    %----------------------G---------------------%
    figure(2);
    phi_str = num2str(phi);
    titlename = sprintf(' NanoParticle : %s \n\\phi = %s',NanoParticle,phi_str);
    
    plot(x,Gu(:,1,i),...
        x,Gu(:,3,i),...
        x,Gu(:,8,i),...
        x,Gu(:,12,i),...
        x,Gu(:,15,i),...
        x,Gs(:,i),'LineWidth',1.2)
    title(titlename);
    xlabel('\eta');
    ylabel('G');
    %legend('t*= 1','2','3','4','5','Steady')
    
    % File name 
    name = sprintf('Data\\%s\\figs\\G\\phi%s',NanoParticle,No);
    
    % Create arrow
    annotation('arrow',[0.16 0.39],...
    [0.18 0.32]);
   
    % Create Text 
    text('String','t*',...
    'Position',[1.73387096774194 0.135479041916168 17.3205080756888],...
    'FontSize',18);
    % Save figure
    print(name,'-djpeg','-r300')
    
    
    %----------------------theta---------------------%
    figure(3);
    phi_str = num2str(phi);
    titlename = sprintf(' NanoParticle : %s \n\\phi = %s',NanoParticle,phi_str);
    
    plot(x,Wu(:,1,i),...
        x,Wu(:,3,i),...
        x,Wu(:,8,i),...
        x,Wu(:,12,i),...
        x,Wu(:,15,i),...
        x,Ws(:,i),'LineWidth',1.2)
    title(titlename);
    xlabel('\eta');
    ylabel('\theta');
    %legend('t*= 1','2','3','4','5','Steady')
    
    % File name 
    name = sprintf('Data\\%s\\figs\\theta\\phi%s',NanoParticle,No);
    
    % Create arrow
    annotation('arrow',[0.14 0.30],...
    [0.17 0.30]);
   
    % Create Text 
    text('String','t*',...
    'Position',[1.16596638655462 0.238911290322581 17.3205080756888],...
    'FontSize',18);
    % Save figure
    print(name,'-djpeg','-r300')   
    
    
    phi = phi + dphi;
 end
 
end
