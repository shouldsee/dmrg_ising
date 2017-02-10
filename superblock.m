function[phi psi]=superblock(T_ml,W)
global oper opts
oddfirst=reshape(reshape(1:8,2,4)',1,[]);
tp=oper(oper(T_ml,shiftdim(T_ml,-4)),permute(W,[5 6 1 2 7 8 3 4]));
ptp=permute(tp,oddfirst);
T_2m=square(ptp);
[Vr,D,Vl]=eigs((T_2m),1,'SA',opts);
% d=diag(D);
id=1;
Vl=Vr;
m=sqrt(max(size(Vl)))/2;
phi=reshape(Vl(:,id)',2*m,1,2*m);
psi=reshape(Vr(:,id),1,2*m,2*m);

end