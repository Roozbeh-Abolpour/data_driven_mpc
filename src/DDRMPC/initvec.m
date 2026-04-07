function v=initvec(Ur,Yr,na,np,nc,nf)
ny=1;nu=1;
nay=na;nau=na;
[a,b]=identifymodel(Ur,Yr,na);
Ass=cell(nc,na);Bss=cell(nc,na);
t=rand(nc,1);t=t/sum(t);
for i=1:nc
    for j=1:nay
        Ass{i,j}=a(j);
    end
    for j=1:nau
        Bss{i,j}=b(j);
    end
end
u=zeros(nf,nu);
y=zeros(nc,nf,ny);
v=t;
for j=1:nay
    for i=1:nc
        v=[v;Ass{i,j}(:)];        
    end
end
for j=1:nau
    for i=1:nc
        v=[v;Bss{i,j}(:)];        
    end
end
v=[v;u(:);y(:)];
end