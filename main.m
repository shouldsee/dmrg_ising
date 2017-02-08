K=0.3;
ss=[-1 1]
mmax=30;

corrf=1;
err=0;
ra=1;

oddfirst=reshape(reshape(1:8,2,4)',1,[]);
[msh1p msh1 msh2p msh2]=ndgrid(ss,ss,ss,ss);
logW=-0.5*K*(msh1.*msh2+msh1p.*msh2p+msh1.*msh1p+msh2.*msh2p);
% W=exp(-0.5*K*(msh1.*msh2+msh1p.*msh2p+msh1.*msh1p+msh2.*msh2p));
% W=exp(logW);
% oper=@times;%     corrf_pre=mean(squeeze(phi).*squeeze(psi),3);
%     corrf=mean(corrf_pre(:));

W=logW;
oper=@plus;
T_ml=W;
figure(2)
fi=imagesc(zeros(2*mmax,2*mmax));
for i=1:1024;

m=size(T_ml,1);
if 4*m<=2*mmax
    % tp=bsxfun(@times,T_ml,shiftdim(W,-2));
    tp=bsxfun(oper,T_ml,shiftdim(W,-2));
    tp=permute(tp,[1 3 2 4 5 6]);
    T_mplr=reshape(tp,2*m,2*m,2,2);
    T_mprr=T_mplr;
    % T_mr=reshape(tp,4*m,4*m,2,2);
else
%% calculate T_2m
%     tp=bsxfun(oper,bsxfun(oper,T_ml,shiftdim(W,-2)),permute(T_ml,[5:8 flip(1:4)]));
%     tp=T_ml+shiftdim(W,-2)+permute(T_ml,[5:8 3 4 1 2]);
%     figure(3)
    tp=T_ml+shiftdim(T_mr,-4)+permute(W,[5 6 1 2 7 8 3 4]);
    T_2m=reshape(permute(tp,oddfirst),sqrt(numel(tp)),[]);
    figure(3)
%     imagesc(T_2m)
%     [Vr,D,Vl]=eigs((T_2m),1);
    [Vr,D,Vl]=eig(T_2m);
    d=diag(D);
%     id=1;
    Vl=Vr;
    [~,id]=max(real([d(1) d(end)]));
%     phi=reshape()
%     Vl=fliplr(Vl);
%     Vr=shiftdim(fliplr(Vr),-1);
%     phi=reshape(Vl(:,id)',2*m,1,2*m);
%     psi=reshape(Vr(:,id),1,2*m,2*m);
    phi=reshape(Vl,2*m,1,2*m,(2*m)^2);
    psi=reshape(Vr,1,2*m,2*m,(2*m)^2);
    
%     corrf_pre=reshape(exp(squeeze(phi)+squeeze(psi)),[m,2,m,2]).*ss.*shiftdim(ss,-2);
    
%     corrf_pre=mean(squeeze(phi).*squeeze(psi),3);
%     corrf=sum(corrf_pre(:));
%     psi=reshape(permute(psi,[1,2,4,3]),1,2*m,[]);
    
    D=rot90(D,2);
    rhoL=squeeze(sum(sum(phi+ psi,3),4));
    rhoR=squeeze(sum(sum(permute(phi,[3 1 2 4])+psi,2),4));
%     imagesc(rhoL)
%     imagesc(rhoR)
    %%
    [T,S,U]=svd(rhoL);
%     [T,S,U]=svds(rhoL,mmax);
    Ol=reshape(T(:,1:mmax)',[mmax,m,2]);
    Ql=shiftdim(reshape(U(:,1:mmax),[m 2 mmax]),-3);
    
    s=diag(S);
%     s(1:5)
%     ra=sum(s(1:mmax))/sum(s);
    
    T_mlo=shiftdim(permute(T_ml,[1 3 2 4]),-1);
    Tst=Ol+T_mlo+Ql;
    T_mpl=Ol + T_mlo +Ql + permute((W),[5,6,1,7,2,8,3,4]);
    T_mpl=squeeze(sum(reshape(T_mpl,mmax,[],mmax,2,2),2));
    %%
    [T,S,U]=svd(rhoR);
%     [T,S,U]=svds(rhoL,mmax);
    Or=reshape(T(:,1:mmax)',[mmax,m,2]);
    Qr=shiftdim(reshape(U(:,1:mmax),[m 2 mmax]),-3);
    
    s=diag(S);
%     s(1:5)
%     ra=sum(s(1:mmax))/sum(s);
    
    T_mro=shiftdim(permute(T_mr,[1 3 2 4]),-1);
%     Tst=Or+T_mro+Qr;
    T_mpr=Or + T_mro +Qr + permute((W),[5,6,1,7,2,8,3,4]);
    T_mpr=squeeze(sum(reshape(T_mpr,mmax,[],mmax,2,2),2));
    
    %%
%     T_mpl=(T_mpl+permute(T_mpl,[2 1 4 3]))/2;
    T_mplr=(T_mpl)-mean(T_mpl(:));
    T_mprr=(T_mpr)-mean(T_mpr(:));
    
end

err=std(T_mplr(:))-std(T_ml(:));
T_ml=T_mplr;
T_mr=T_mprr;
fprintf('error is %.4d ,corrf is %.4f, ratio is %.4f \n ',err,corrf,ra);    
vsa=reshape(permute(T_ml,[1 3 2 4]),numel(T_ml)^0.5,[]);
siz=size(vsa);
vs=padarray(vsa,[2*mmax 2*mmax]-siz,0,'post');
% set(fi,'CData',real(vs));
% drawnow
% pause(0.5)
end% T_ml=
%%
imagesc(reshape(T_ml,numel(T_ml)^0.5,[]))
%%
imagesc(reshape(T_mpl,numel(T_ml)^0.5,[]))
