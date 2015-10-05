clear 
clc

p = 0 ;
dp = 0.01;
MAXp = 0.21;
NanoParticle = 'Cu';

F1=zeros(100,21);
G1=zeros(100,21);
W1=zeros(100,21);
m=1;
totaltime = 0;
while (p <= MAXp)
    No = strrep(num2str(p), '.', '_');
    name = sprintf('Data\\%s\\phi%s.mat',NanoParticle,No);
    load(name)
    time
    totaltime=totaltime + time;
    F1(:,m)=F(:,2);
    G1(:,m)=G(:,2);
    W1(:,m)=W(:,2);
    m = m+1;
    p = p + dp;
end
figure(1)
plot(x,F1(:,1),x,F1(:,7),x,F1(:,14)...
    ,x,F1(:,21));
xlabel('n')
ylabel('F''');
legend('0','0.06','0.13','0.2')

figure(2)
plot(x,G1(:,1),x,G1(:,7),x,G1(:,14)...
    ,x,G1(:,21));
xlabel('n')
ylabel('G''');
legend('0','0.06','0.13','0.2')

figure(3)
plot(x,W1(:,1),x,W1(:,3),x,W1(:,7),x,W1(:,10),x,W1(:,14)...
    ,x,W1(:,21));
xlabel('n')
ylabel('W''');
legend('0','0.02','0.06','0.09','0.13','0.2')
totaltime
totaltime/60