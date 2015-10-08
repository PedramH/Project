%Puts all important data in a single mat file
clear 
clc

p = 0 ;
dp = 0.01;
MAXp = 0.31;
NanoParticle = 'Cu';

FileNo = 31;
n=100;
j=50;

Fs=zeros(n,FileNo);        %Steady state data
Gs=zeros(n,FileNo);
Ws=zeros(n,FileNo);

Fu=zeros(n,j,FileNo);        %unsteady data
Gu=zeros(n,j,FileNo);
Wu=zeros(n,j,FileNo);

Kn=zeros(FileNo);      %Knf


m=1;


while (p <= MAXp)
    No = strrep(num2str(p), '.', '_');
    name = sprintf('Data\\%s\\phi%s.mat',NanoParticle,No);
    load(name)
    Fs(:,m)=F(:,2);
    Gs(:,m)=G(:,2);
    Ws(:,m)=W(:,2);
    
    Fu(:,:,m)=Fr(:,:);
    Gu(:,:,m)=Gr(:,:);
    Wu(:,:,m)=Wr(:,:);
    
    Kn(m) = Knf;
    
    m = m+1;
    p = p + dp;
end
clear Knf
Knf = zeros(FileNo);
Knf(:) = Kn(:);
name = sprintf('Data\\%s\\AllData.mat',NanoParticle);
save(name,'Fs','Gs','Ws','Fu','Gu','Wu','MAXphi','dphi',...
    'NanoParticle','x','C','Kf','Knf','FileNo','n')