function y=ddnet(U,Y,Wu,Wy,b,sig)
y=sig(Wu*U+Wy*Y+b);
end