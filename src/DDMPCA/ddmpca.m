function [uo,ef]=ddmpca(u,y,yref,bounds,orders,k)
global U Y ab;
umin=bounds.umin;umax=bounds.umax;
ymin=bounds.ymin;ymax=bounds.ymax;
dumax=bounds.dumax;
L=orders.L;n=orders.n;N=orders.N;
Nc=N;
if k==0
    U=zeros(L+n,Nc);
    Y=zeros(L+n,Nc);
    uo=0.1*(umin+rand*(umax-umin));
    ab=zeros(Nc,1);
    ef=0;
    return
end
U1=U(:,2:end);U2=[U(2:end,end);u];
U=[U1 U2];
Y1=Y(:,2:end);Y2=[Y(2:end,end);y];
Y=[Y1 Y2];
if k<=L+N+n
    uo=0.1*(umin+rand*(umax-umin));
    ef=0;
    return
end
e=ones(L,1);et=zeros(1,L);et(1)=1;
It=eye(L)-diag(ones(1,L-1),1);It(end,:)=[];

Umin=e*umin;Umax=e*umax;
Ymin=e*ymin;Ymax=e*ymax;
Dumax=ones(L-1,1)*dumax;
Yd=e*yref;

Ue=U(n+1:end,:);
Ye=Y(n+1:end,:);
H=2*(Ye'*Ye);mh=min(eig(H));
if mh<0
    H=H+1e-4*eye(Nc);
end
f=-2*Ye'*Yd;

Aeq=[U(1:n,:);Y(1:n,:)];
beq=[U(end-n+1:end,end);Y(end-n+1:end,end)];

A=[-Ue;Ue;-Ye;Ye;It*Ue;-It*Ue;et*Ue;-et*Ue];
b=[-Umin;Umax;-Ymin;Ymax;Dumax;Dumax;u+dumax;-u+dumax];

op=optimset('disp','none');
try
[a,~,ef]=quadprog(H,f,A,b,Aeq,beq,[],[],ab,op);
catch
    uo=0.1*(umin+rand*(umax-umin));
    ef=0;
    return
end
if ~isempty(a)&&ef>=0
    ab=a;
    uc=Ue(1,:)*a;    
    ef=1;
else
    uc=0.1*(umin+rand*(umax-umin));
    ef=0;
end
uo=uc;
end