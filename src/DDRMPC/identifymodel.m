function [a,b]=identifymodel(Ur,Yr,na)
nr=length(Ur);nv=2*na;
it=1;inds=1:nv;
IA=reshape(inds(it:it+na-1),na,1);it=it+na;
IB=reshape(inds(it:it+na-1),na,1);
A=[];b=[];
for i=1:nr-na
    a=zeros(1,nv);
    for j=1:na
        a(IA(j))=Yr(i+j);
        a(IB(j))=Ur(i+j-1);
    end
    A=[A;a];b=[b;Yr(i)];
end
th=A\b;
a=th(1:end/2);b=th(end/2+1:end);
end