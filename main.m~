K=-0.1;
ss=[-1 1]
mmax=16;

oddfirst=reshape(reshape(1:8,2,4)',1,[]);
[msh1p msh1 msh2p msh2]=ndgrid(ss,ss,ss,ss);
logW=-0.5*K*(msh1.*msh2+msh1p.*msh2p+msh1.*msh1p+msh2.*msh2p);
% W=exp(-0.5*K*(msh1.*msh2+msh1p.*msh2p+msh1.*msh1p+msh2.*msh2p));
% W=exp(logW);
% oper=@times;
W=logW;
oper=@plus;
T_ml=W;
figure(2)
fi=imagesc(zeros(2*mmax,2*mmax));
for i=1:50;

m=size(T_ml,1);
if 4*m<=2*mmax
    % tp=bsxfun(@times,T_ml,shiftdim(W,-2));
    tp=bsxfun(oper,T_ml,shiftdim(W,-2));
    tp=permute(tp,[1 3 2 4 5 6]);
    T_mplr=reshape(tp,2*m,2*m,2,2);
    % T_mr=reshape(tp,4*m,4*m,2,2);
else
%% calculate T_2m
%     tp=bsxfun(oper,bsxfun(oper,T_ml,shiftdim(W,-2)),permute(T_ml,[5:8 flip(1:4)]));
%     tp=T_ml+shiftdim(W,-2)+permute(T_ml,[5:8 3 4 1 2]);
%     figure(3)
    tp=T_ml+shiftdim(T_ml,-4)+permute(W,[5 6 1 2 7 8 3 4]);
    tpr=reshape(permute(tp,oddfirst),sqrt(numel(tp)),[]);
    
    [Vr,D,Vl]=eig((tpr));
    d=diag(D);
    [~,id]=max(real([d(1) d(end)]));
%     phi=reshape()
%     Vl=fliplr(Vl);
%     Vr=shiftdim(fliplr(Vr),-1);
    phi=reshape(Vl(:,id)',2*m,1,2*m);
    psi=reshape(Vr(:,id),1,2*m,2*m);
    corrf_pre=reshape(squeeze(phi).*squeeze(psi),[m,2,m,2]).*ss.*shiftdim(ss,-2);
    corrf=sum(corrf_pre(:));
%     psi=reshape(permute(psi,[1,2,4,3]),1,2*m,[]);
    
    D=rot90(D,2);
    rhoL=squeeze(sum(phi+ psi,3));
    
%     imagesc(rhoL)
    %%
    [T,S,U]=svd(rhoL);
    Ol=reshape(T(:,1:mmax)',[mmax,m,2]);
    Ql=shiftdim(reshape(U(:,1:mmax),[m 2 mmax]),-3);
    s=diag(S);
    ra=sum(s(1:mmax))/sum(s);
    
    T_mlo=shiftdim(permute(T_ml,[1 3 2 4]),-1);
    Tst=Ol+T_mlo+Ql;
    T_mpl=Ol + T_mlo +Ql + permute((W),[5,6,1,7,2,8,3,4]);
    T_mpl=squeeze(sum(reshape(T_mpl,m,[],m,2,2),2));
    %%
    
    T_mplr=(T_mpl)-mean(T_mpl(:));
    
end

err=std(T_mplr(:))-std(T_ml(:));
T_ml=T_mplr;
fprintf('error is %.4d ,corrf is %.4f, ratio is %.4f \n ',err,corrf,ra);    
vsa=reshape(permute(T_ml,[1 3 2 4]),numel(T_ml)^0.5,[]);
siz=size(vsa);
vs=padarray(vsa,[2*mmax 2*mmax]-siz,0,'post');
set(fi,'CData',vs);
drawnow
pause(0.5)
end% T_ml=
%%
imagesc(reshape(T_ml,numel(T_ml)^0.5,[]))
%%
imagesc(reshape(T_mpl,numel(T_ml)^0.5,[]))
