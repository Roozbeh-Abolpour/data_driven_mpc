function [uc,ef]=kfddmpc(u,y,yref,bounds,orders,k)
global U Y P th Uo yh;
umin=bounds.umin;umax=bounds.umax;
ymin=bounds.ymin;ymax=bounds.ymax;
dumax=bounds.dumax;

Nu=orders.Nu;Ny=orders.Ny;
p=orders.p;q=orders.q;

if k==0      
    P=1e-2*eye(p+q);th=zeros(p+q,1);
    U=zeros(q,1);Y=zeros(p,1);
    Uo=zeros(Nu,1);    
end
U=[U(2:end);u];

H=[Y(p:-1:1);U(q:-1:1)];
L=(P*H)/(1+H'*P*H);
P=P-L*H'*P;
th=th+L*(y-H'*th);
yh=0;
for i=1:p        
    yh=yh+th(i)*Y(p-i+1,:)';
end
for i=1:q  
    yh=yh+th(i+p)*U(q-i+1,:)';
end
Y=[Y(2:end);y];

if k>2*(p+q)
G=zeros(Ny,Nu);g=zeros(Ny,1);
for i=1:Ny
    Ge=zeros(1,size(G,2));ge=0;
    for j=1:p        
        aj=th(j);
        if j>=i
            ge=ge+aj*Y(end-j+i,:)';
        else
            Gw=G(i-j,:);gw=g(i-j,:);
            Ge=Ge+aj*Gw;ge=ge+aj*gw;
        end
    end
    for j=1:q        
        bj=th(p+j);
        if j>i
            ge=ge+bj*U(end-j+i+1,:)';
        elseif i-j+1<=Nu
            Ge(:,i-j+1)=Ge(:,i-j+1)+bj;
        end
    end
    G(i,:)=Ge;
    g(i)=ge;
end
z1=ones(Nu,1);z2=ones(Ny,1);
Umin=z1*umin;Umax=z1*umax;
Ymin=z2*ymin;Ymax=z2*ymax;
Dumax=z1(1:end-1)*dumax;

et=zeros(1,Nu);et(1)=1;
It=eye(Nu)-diag(ones(1,Nu-1),1);It(end,:)=[];

A=[-G;G;It;-It;et;-et];
b=[g-Ymin;-g+Ymax;Dumax;Dumax;u+dumax;-u+dumax];

Yref=z2*yref;gb=g-Yref;
opt=optimset('disp','none');
[Ue,~,ef]=quadprog(2*(G'*G),2*G'*gb,A,b,[],[],Umin,Umax,Uo,opt);
if ~isempty(Ue)&&ef>0&&k>max(Ny,Nu)
    Uo=Ue;
    ef=1;
    uc=Uo(1);
else
    ef=0;
    uc=0.1*(umin+rand*(umax-umin));
end
else
    ef=0;
    uc=0.1*(umin+rand*(umax-umin));

end
end