function u=ddrmpc(u,y,yref,bounds,orders,k)
global Ur Yr v g;
yd=yref;nu=1;ny=1;
umin=bounds.umin;umax=bounds.umax;
ymin=bounds.ymin;ymax=bounds.ymax;
udmax=bounds.dumax;

na=orders.na;nf=orders.nf;
nc=orders.nc;np=orders.np;

if k==0
    Ur=zeros(nu,np);Yr=zeros(ny,np);
end
ub=u;
Ur(:,2:end)=Ur(:,1:end-1);Ur(:,1)=u;
Yr(:,2:end)=Yr(:,1:end-1);Yr(:,1)=y;

if k==np
    v=initvec(Ur,Yr,na,np,nc,nf);
    g=10;
    v=[v;g];
end
if k>np
    nay=na;nau=na;nr=np-max(na);
    nv=1+nc+nay*nc*ny*ny+nau*nc*ny*nu+nf*nu+nc*nf*ny;
    it=1;inds=1:nv;
    It=reshape(inds(it:it+nc-1),nc,1);it=it+nc;
    IA=reshape(inds(it:it+nc*nay*ny*ny-1),nc,nay,ny,ny);it=it+nc*nay*ny*ny;
    IB=reshape(inds(it:it+nc*nau*ny*nu-1),nc,nau,ny,nu);it=it+nc*nau*ny*nu;
    Iu=reshape(inds(it:it+nf*nu-1),nf,nu);it=it+nf*nu;
    Iy=reshape(inds(it:it+nc*nf*ny-1),nc,nf,ny);it=it+nc*nf*ny;
    Ig=it;
    Qs={};
    for i=1:nr

        Q=zeros(nv+1);
        for s=1:nay
            for l=1:nc

                Q(IA(l,s,1,1),It(l))=0.5*Yr(1,i+s);
                Q(It(l),IA(l,s,1,1))=0.5*Yr(1,i+s);

            end
        end
        for s=0:nau-1
            for l=1:nc

                Q(IB(l,s+1,1,1),It(l))=0.5*Ur(1,i+s);
                Q(It(l),IB(l,s+1,1,1))=0.5*Ur(1,i+s);

            end
        end
        Qb=Q;
        Q(end,end)=-Yr(1,i);
        Q(end,Ig)=-0.5;Q(Ig,end)=-0.5;
        Qs{end+1}=Q;
        Q=-Qb;
        Q(end,end)=+Yr(1,i);
        Q(end,Ig)=-0.5;Q(Ig,end)=-0.5;
        Qs{end+1}=Q;

    end
    for l=1:nc
        for i=1:nf

            Q=zeros(nv+1);
            for s=1:nay

                if i<=s
                    Q(IA(l,s,1,1),end)=0.5*Yr(1,s-i+1);
                    Q(end,IA(l,s,1,1))=0.5*Yr(1,s-i+1);
                else
                    Q(IA(l,s,1,1),Iy(l,i-s,1))=0.5;
                    Q(Iy(l,i-s,1),IA(l,s,1,1))=0.5;
                end

            end
            for s=0:nau-1
                if i<=s
                    Q(IB(l,s+1,1,1),end)=0.5*Ur(1,s-i+1);
                    Q(end,IB(l,s+1,1,1))=0.5*Ur(1,s-i+1);
                else
                    Q(IB(l,s+1,1,1),Iu(i-s,1))=0.5;
                    Q(Iu(i-s,1),IB(l,s+1,1,1))=0.5;
                end

            end
            Q(end,Iy(l,i,1))=-0.5;Q(Iy(l,i,1),end)=-0.5;
            Qb=Q;
            Q(end,end)=-1e-2;
            Qs{end+1}=Q;
            Q=-Qb;
            Q(end,end)=-1e-2;
            Qs{end+1}=Q;
        end

    end
    A=[];b=[];
    for l=1:nc
        for i=1:nf
            for j=1:ny
                a=zeros(1,nv);
                a(Iy(l,i,j))=1;a(Ig)=-1;
                A=[A;a];b=[b;yd(j)];
                a=zeros(1,nv);
                a(Iy(l,i,j))=-1;a(Ig)=-1;
                A=[A;a];b=[b;-yd(j)];
            end
        end
    end
    for i=1:nf
        for j=1:nu
            a=zeros(1,nv);
            a(Iu(i,j))=1;
            A=[A;a];b=[b;umax(j)];
            a=zeros(1,nv);
            a(Iu(i,j))=-1;
            A=[A;a];b=[b;-umin(j)];
        end
    end
    for l=1:nc
        for i=1:nf
            for j=1:ny
                a=zeros(1,nv);
                a(Iy(l,i,j))=1;
                A=[A;a];b=[b;ymax(j)];
                a=zeros(1,nv);
                a(Iy(l,i,j))=-1;
                A=[A;a];b=[b;-ymin(j)];
            end
        end
    end
    for j=1:nu
        a=zeros(1,nv);
        a(Iu(1,j))=1;
        A=[A;a];b=[b;ub(j)+udmax(j)];
        a=zeros(1,nv);
        a(Iu(1,j))=-1;
        A=[A;a];b=[b;-ub(j)+udmax(j)];
    end
    for i=2:nf
        for j=1:nu
            a=zeros(1,nv);
            a(Iu(i,j))=1;a(Iu(i-1,j))=-1;
            A=[A;a];b=[b;udmax(j)];
            a=zeros(1,nv);
            a(Iu(i,j))=-1;a(Iu(i-1,j))=1;
            A=[A;a];b=[b;udmax(j)];
        end
    end
    for l=1:nc
        a=zeros(1,nv);
        a(It(l))=-1;
        A=[A;a];b=[b;0];
    end
    a=zeros(1,nv);
    a(It)=1;
    A=[A;a];b=[b;1+1e-2];
    a=zeros(1,nv);
    a(It)=-1;
    A=[A;a];b=[b;-1+1e-2];


    % for i=1:length(b)
    %     Q=zeros(nv+1);
    %     Q(1:end-1,end)=0.5*A(i,:)';
    %     Q(end,1:end-1)=0.5*A(i,:);
    %     Q(end,end)=-b(i);
    %     Qs{end+1}=Q;
    % end
    % A=[];b=[];
    [Qhs,Ah,kapa]=qcfp2fieos(Qs,A,b,1e-2);
    vh=fmap(v,kapa,Qhs);
    vh=ieosstep(Qhs,Ah,vh);
    v=rmap(vh,kapa);
    if vh'*vh<1-1e-2||(abs(y-mean(v(Iy(:))))>0.9&&g<=0.01)
        u=0.01*(umin+rand*(umax-umin));
        v=initvec(Ur,Yr,na,np,nc,nf);
        g=1;
        v=[v;g];
    else
        u=v(Iu(1,1));
        g=0.9*g;
    end    
    v(end)=g;
else
    u=0.1*(umin+rand*(umax-umin));
end
end