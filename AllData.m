%Puts all important data in a single mat file
%First index of all martices is the NanoParticle index
clear 
clc

p = 0 ;
dp = 0.01;
MAXp = 0.21;
NanoParticles = {'Cu','Al2O3','TiO2'};


FileNo = 21;
n=100;
j=350;

Fs=zeros(3,n,FileNo);        %Steady state data
Gs=zeros(3,n,FileNo);
Ws=zeros(3,n,FileNo);

Fu=zeros(3,n,j,FileNo);        %unsteady data
Gu=zeros(3,n,j,FileNo);
Wu=zeros(3,n,j,FileNo);

Kn=zeros(3,FileNo);      %Knf


for in=1:3
    
    NanoParticle = NanoParticles{in}
    m=1;
    p = 0 ;
    while (p <= MAXp)
        No = strrep(num2str(p), '.', '_');
        name = sprintf('Data\\%s\\phi%s.mat',NanoParticle,No);
        load(name)
        Fs(in,:,m)=F(:,2);
        Gs(in,:,m)=G(:,2);
        Ws(in,:,m)=W(:,2);
    
        
        
        Fu(in,:,1:size(Fr,2),m)=Fr(:,:);
        Gu(in,:,1:size(Gr,2),m)=Gr(:,:);
        Wu(in,:,1:size(Wr,2),m)=Wr(:,:);
    
        Kn(in,m) = Knf;
    
        m = m+1;
        p = p + dp;
    end
    
end
    clear Knf
    Knf = zeros(3,FileNo);
    Knf(:,:) = Kn(:,:);
    name = sprintf('Data\\AllData.mat');
    save(name,'Fs','Gs','Ws','Fu','Gu','Wu','MAXphi','dphi',...
        'NanoParticles','x','C','Kf','Knf','FileNo','n','detta','dt')