function[corrf]=NNcorr(phi,psi,ss)
global oper m
p=oper(squeeze(phi),squeeze(psi));        
corrf_pre=reshape(p,[m,2,m,2]).*(ss.*shiftdim(ss,-2));....*shiftdim(ss,-2);
corrf=sum(corrf_pre(:));

end