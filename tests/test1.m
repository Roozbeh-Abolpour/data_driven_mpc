close all
clear
clc

f=@(x,u)[0.1*x(1)+x(1)*x(2);0.5*x(2)-sin(x(1))+cos(0.1*x(2))*u];
g=@(x,u) x(1)+x(2);
N=100;
bounds.umin=-2;bounds.umax=2;
bounds.ymin=-10;bounds.ymax=10;
bounds.dumax=0.1;

yrefs=ones(1,N);
umins=bounds.umin*ones(1,N);umaxs=bounds.umax*ones(1,N);
ymins=bounds.ymin*ones(1,N);ymaxs=bounds.ymax*ones(1,N);

ordersp.na=3;ordersp.nf=2;ordersp.nc=1;ordersp.np=4;
ordersf.Nu=10;ordersf.Ny=10;ordersf.p=2;ordersf.q=2;
ordersa.L=10;ordersa.n=2;ordersa.N=40;
ordersn.nd=20;ordersn.nf=5;

x=[1;-1];
u=0;y=0;
ups=zeros(1,N);
yps=zeros(1,N);
for k=0:N-1
    u=ddrmpc(u,y,yrefs(k+1),bounds,ordersp,k);
    x=f(x,u);
    y=g(x,u);
    ups(k+1)=u;
    yps(k+1)=y;
end

x=[1;-1];
u=0;y=0;
ufs=zeros(1,N);
yfs=zeros(1,N);
for k=0:N-1
    u=kfddmpc(u,y,yrefs(k+1),bounds,ordersf,k);
    x=f(x,u);
    y=g(x,u);
    ufs(k+1)=u;
    yfs(k+1)=y;
end

x=[1;-1];
u=0;y=0;
uas=zeros(1,N);
yas=zeros(1,N);
for k=0:N-1
    u=ddmpca(u,y,yrefs(k+1),bounds,ordersa,k);
    x=f(x,u);
    y=g(x,u);
    uas(k+1)=u;
    yas(k+1)=y;
end

x=[1;-1];
u=0;y=0;
uns=zeros(1,N);
yns=zeros(1,N);
for k=0:N-1
    u=nnddmpc(u,y,yrefs(k+1),bounds,ordersn,k);
    x=f(x,u);
    y=g(x,u);
    uns(k+1)=u;
    yns(k+1)=y;
end

figure
hold on
grid on
box on
plot(1:N,yps,'LineWidth',1.5)
plot(1:N,yfs,'LineWidth',1.5)
plot(1:N,yns,'LineWidth',1.5)
plot(1:N,yas,'LineWidth',1.5)
plot(1:N,yrefs,'g--','LineWidth',1.5)
plot(1:N,ymins,'r:','LineWidth',1.2)
plot(1:N,ymaxs,'r:','LineWidth',1.2)
xlabel('time step (k)')
ylabel('system output (y)')
title('Output comparison')
legend('DD-RMPC','KFDDMPC','NNDDMPC','DDMPCA','reference','bounds','Location','best')
set(gca,'fontname','timesnewroman')
set(gca,'fontsize',13)

figure
hold on
grid on
box on
plot(1:N,ups,'LineWidth',1.5)
plot(1:N,ufs,'LineWidth',1.5)
plot(1:N,uns,'LineWidth',1.5)
plot(1:N,uas,'LineWidth',1.5)
plot(1:N,umins,'r:','LineWidth',1.2)
plot(1:N,umaxs,'r:','LineWidth',1.2)
xlabel('time step (k)')
ylabel('system input (u)')
title('Input comparison')
legend('DD-RMPC','KFDDMPC','NNDDMPC','DDMPCA','bounds','Location','best')
set(gca,'fontname','timesnewroman')
set(gca,'fontsize',13)

figure
hold on
grid on
box on
plot(1:N,yps,'LineWidth',1.5)
plot(1:N,yrefs,'g--','LineWidth',1.5)
plot(1:N,ymins,'r:','LineWidth',1.2)
plot(1:N,ymaxs,'r:','LineWidth',1.2)
xlabel('time step (k)')
ylabel('system output (y)')
title('DD-RMPC output')
legend('DD-RMPC','reference','bounds','Location','best')
set(gca,'fontname','timesnewroman')
set(gca,'fontsize',13)

figure
hold on
grid on
box on
plot(1:N,ups,'LineWidth',1.5)
plot(1:N,umins,'r:','LineWidth',1.2)
plot(1:N,umaxs,'r:','LineWidth',1.2)
xlabel('time step (k)')
ylabel('system input (u)')
title('DD-RMPC input')
legend('DD-RMPC','bounds','Location','best')
set(gca,'fontname','timesnewroman')
set(gca,'fontsize',13)

figure
hold on
grid on
box on
plot(1:N,yfs,'LineWidth',1.5)
plot(1:N,yrefs,'g--','LineWidth',1.5)
plot(1:N,ymins,'r:','LineWidth',1.2)
plot(1:N,ymaxs,'r:','LineWidth',1.2)
xlabel('time step (k)')
ylabel('system output (y)')
title('KFDDMPC output')
legend('KFDDMPC','reference','bounds','Location','best')
set(gca,'fontname','timesnewroman')
set(gca,'fontsize',13)

figure
hold on
grid on
box on
plot(1:N,ufs,'LineWidth',1.5)
plot(1:N,umins,'r:','LineWidth',1.2)
plot(1:N,umaxs,'r:','LineWidth',1.2)
xlabel('time step (k)')
ylabel('system input (u)')
title('KFDDMPC input')
legend('KFDDMPC','bounds','Location','best')
set(gca,'fontname','timesnewroman')
set(gca,'fontsize',13)
figure
hold on
grid on
box on
plot(1:N,yns,'LineWidth',1.5)
plot(1:N,yrefs,'g--','LineWidth',1.5)
plot(1:N,ymins,'r:','LineWidth',1.2)
plot(1:N,ymaxs,'r:','LineWidth',1.2)
xlabel('time step (k)')
ylabel('system output (y)')
title('NNDDMPC output')
legend('NNDDMPC','reference','bounds','Location','best')
set(gca,'fontname','timesnewroman')
set(gca,'fontsize',13)

figure
hold on
grid on
box on
plot(1:N,uns,'LineWidth',1.5)
plot(1:N,umins,'r:','LineWidth',1.2)
plot(1:N,umaxs,'r:','LineWidth',1.2)
xlabel('time step (k)')
ylabel('system input (u)')
title('NNDDMPC input')
legend('NNDDMPC','bounds','Location','best')
set(gca,'fontname','timesnewroman')
set(gca,'fontsize',13)
figure
hold on
grid on
box on
plot(1:N,yas,'LineWidth',1.5)
plot(1:N,yrefs,'g--','LineWidth',1.5)
plot(1:N,ymins,'r:','LineWidth',1.2)
plot(1:N,ymaxs,'r:','LineWidth',1.2)
xlabel('time step (k)')
ylabel('system output (y)')
title('DDMPCA output')
legend('DDMPCA','reference','bounds','Location','best')
set(gca,'fontname','timesnewroman')
set(gca,'fontsize',13)
figure
hold on
grid on
box on
plot(1:N,uas,'LineWidth',1.5)
plot(1:N,umins,'r:','LineWidth',1.2)
plot(1:N,umaxs,'r:','LineWidth',1.2)
xlabel('time step (k)')
ylabel('system input (u)')
title('DDMPCA input')
legend('DDMPCA','bounds','Location','best')
set(gca,'fontname','timesnewroman')
set(gca,'fontsize',13)

