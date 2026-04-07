function [Qhs,Ah,kapa]=qcfp2fieos(Qs,A,b,kapa0)
m=length(Qs);
for i=1:m
    Qs{i}(end,end)=Qs{i}(end,end)-1;
end
kapa=1;n=length(Qs{1})-1;
for i=1:m
    w=min(eig(Qs{i}));    
    if w<0
        kapa=max(kapa,-w);
    end
end
kapa=kapa+kapa0;
q1=1/sqrt(kapa);q2=1/sqrt(kapa-1);
q3=sqrt(kapa-1)/sqrt(kapa);
D=diag([q1*ones(1,n) q2]);
Qhs=Qs;
for i=1:length(Qs)
    Qhs{i}=D'*(Qs{i}+kapa*eye(n+1))*D;
end
if ~isempty(A)
    Ah=[q3*A -b];    
else
    Ah=[];
end
end