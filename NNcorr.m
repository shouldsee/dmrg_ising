function[corrf,mag]=NNcorr(phi,psi,ss)
global oper m
p=reshape(oper(squeeze(phi),squeeze(psi)),[m 2 m 2]);        
corrf_pre2=p.*(ss.*shiftdim(ss,-2));....*shiftdim(ss,-2);
% mag_pre=p.*(0.5*(ss+shiftdim(ss,-2))).^2;
ssj=shiftdim(ss,-2);
% msite=(ss.^2+ssj.^2-(ss+ssj).^2)/2;
msite=ss+ssj;
mag_pre=p.*msite;
mag=sum(mag_pre(:));
% mag_pre=reshape(p,[m,2,m,2]).*(ss+shiftdim(ss,-2));
% mag=sum(mag_pre(:));
% corrf_pre1a=reshape(p,[m,2,m,2]).*(ss);....*shiftdim(ss,-2);
% corrf_pre1b=reshape(p,[m,2,m,2]).*shiftdim(ss,-2);....*shiftdim(ss,-2);
% corrf=sum(corrf_pre2(:));
% corrf=sum(corrf_pre2(:))-sum(corrf_pre1a(:))*sum(corrf_pre1b(:));
corrf=sum(corrf_pre2(:));

end