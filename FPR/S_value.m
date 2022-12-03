function [S,p,r,F]=S_value(M,L,E,S0)
Mask= ForegroundMask(E,M,L);Mask=double(Mask);
S0(S0~=0)=1;S0=double(S0);
c=rdivide(Mask,S0);TP=sum(c(:)==1);
cc=S0-Mask;FN=sum(cc(:)==1);
ccc=Mask-S0;FP=sum(ccc(:)==1);
S=Mask;p=TP/(TP+FP);r=TP/(TP+FN);F=2*p*r/(p+r);
end