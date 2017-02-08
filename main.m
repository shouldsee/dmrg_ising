K=0.44;
ss=[-1 1]
mmax=30;

corrf=1;
err=1;
ra=1;

    opts.disp = 0;   % disable diagnostic information in eigs
    opts.issym = 1;
    opts.real = 1;
%     [Psi Energy] = eigs(H_super,1,'SA', opts);


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
corrs=[];
% for i=1:1024;
% Ts=0.5:0.1:3;
Ks=0.4:0.01:0.5;
% Ts=1./Ks;
Ks=0.804;
for K=Ks;
T=1/K;
err=1;
logW=-0.5*K*(msh1.*msh2+msh1p.*msh2p+msh1.*msh1p+msh2.*msh2p);
W=logW;
T_ml=W;
    while abs(err)>1e-15 | 2*m<mmax;
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
        tp=T_ml+shiftdim(T_ml,-4)+permute(W,[5 6 1 2 7 8 3 4]);
    %     tp=T_ml+shiftdim(T_mr,-4)+permute(W,[5 6 1 2 7 8 3 4]);
        ptp=permute(tp,oddfirst);
    %     aptp=square(ptp); 
    %     sum(sum(aptp-aptp'))
    %     sptp=0.5*(ptp+permute(ptp,[5 6 7 8 1 2 3 4]));
    %     T_2m=reshape(sptp,sqrt(numel(tp)),[]);
    T_2m=square(ptp);
    %     figure(3)
    %     imagesc(T_2m)
        [Vr,D,Vl]=eigs((T_2m),1,'LA',opts);
    
    %     [Vr,D,Vl]=eig(T_2m);
        d=diag(D);
    %     id=1;
        Vl=Vr;
        [~,id]=min(real([d(1) d(end)]));
    %     phi=reshape()
    %     Vl=fliplr(Vl);
    %     Vr=shiftdim(fliplr(Vr),-1);
        phi=reshape(Vl(:,id)',2*m,1,2*m);
        psi=reshape(Vr(:,id),1,2*m,2*m);
        p=squeeze(phi).*squeeze(psi);
        corrf_pre=reshape(p,[m,2,m,2]).*(ss.*shiftdim(ss,-2));....*shiftdim(ss,-2);
%         corrf_pre=mean(squeeze(phi).*squeeze(psi),3);
%         corrf=2*K*sum(corrf_pre(:));
        corrf=sum(corrf_pre(:));
    %     psi=reshape(permute(psi,[1,2,4,3]),1,2*m,[]);

        D=rot90(D,2);
        rhoL=squeeze(sum(sum(phi+ psi,3),4));
    %     rhoR=squeeze(sum(sum(permute(phi,[3 1 2 4])+psi,2),4));
    %     imagesc(rhoL)
    %     imagesc(rhoR)
        %%

%         [V S] = eigs(rhoL,mmax,'LA');
%         ra=sum(diag(S))/trace(abs(rhoL));
% %         [D Index] = sort(diag(D), 'descend');  % sort eigenvalues descending
% %         V = V(:,Index);   
%         Ol=reshape(V',[mmax m 2]);
%         Ql=shiftdim(reshape(V,[m 2 mmax]),-3);

    %     [T,S,U]=svd(rhoL);
%         [T,S,U]=svds(rhoL,mmax);
%         Ol=reshape(T(:,1:mmax)',[mmax,m,2]);
%         Ql=shiftdim(reshape(U(:,1:mmax),[m 2 mmax]),-3);
    %     
    %     s=diag(S);
    % %     s(1:5)
    % %     ra=sum(s(1:mmax))/sum(s);
    %     
%         T_mlo=shiftdim(permute(T_ml,[1 3 2 4]),-1);
    %     a=T_ml+shiftdim(W,-2);
    %     b=shiftdim(reshape(permute(a,[1 3 2 4 5 6]),2*m,2*m,1,2,2),-1);
    %     T_mpl=squeeze(sum(sum(V'+b+shiftdim(V,-2),2),3));
    %     T_mpl=0.5*(T_mpl+permute(T_mpl,[2 1 4 3]));
    %     
%     Tst=sum(reshape(Ol+T_mlo+Ql,mmax,[],mmax),2);

if strcmp(opt.decomp,'svd');
        %%rho=U*D*D*U';
%         [U,S,T]=svds(squeeze(psi),mmax);
%         imagesc(T'*S.^2*T);
%         imagesc(rhoL)
        [T,S,U]=svds(rhoL,mmax);
        Ol=reshape(T(:,1:mmax)',[mmax,m,2]);
        Ql=shiftdim(reshape(U(:,1:mmax),[m 2 mmax]),-3);
        ra=sum(diag(S))/trace(abs(rhoL));
elseif strcmp(opt.decomp,'eig');
        [V S] = eigs(rhoL,mmax,'LA');
        pV=pinv(V);
        ra=sum(diag(S))/trace(abs(rhoL));
%         [D Index] = sort(diag(D), 'descend');  % sort eigenvalues descending
%         V = V(:,Index);   
        Ol=reshape(pV,[mmax m 2]);
        Ql=shiftdim(reshape(V,[m 2 mmax]),-3);
end
    T_mlo=shiftdim(permute(T_ml,[1 3 2 4]),-1);
    T_mpl=Ol + T_mlo +Ql + permute((W),[5,6,1,7,2,8,3,4]);
    T_mpl=squeeze(sum(reshape(T_mpl,mmax,[],mmax,2,2),2));
%     T_mpl=reshape(permute(T_mpl,[1 6 7 8 2 3 4 5]),mmax,mmax,2,2,[] );
%     T_mpl=squeeze(sum(T_mpl,5)); 
%     T_mplv=permute(T_mpl,[1 3 2 4]);
        T_mpl=0.5*(T_mpl+permute(T_mpl,[2 1 4 3]));
        T_mplr=T_mpl-mean(T_mpl(:));
%         T_mplr=T_mpl;
        
    %     %     T_mpl=squeeze(sum(reshape(T_mpl,mmax,[],mmax,2,2),2));
        %%
%         [T,S,U]=svd(rhoR);
    % %     [T,S,U]=svds(rhoL,mmax);
    %     Or=reshape(T(:,1:mmax)',[mmax,m,2]);
    %     Qr=shiftdim(reshape(U(:,1:mmax),[m 2 mmax]),-3);
    %     
    %     s=diag(S);
    % %     s(1:5)
    % %     ra=sum(s(1:mmax))/sum(s);
    %     
    %     T_mro=shiftdim(permute(T_mr,[1 3 2 4]),-1);
    % %     Tst=Or+T_mro+Qr;
    %     T_mpr=Or + T_mro +Qr + permute((W),[5,6,1,7,2,8,3,4]);
    %     T_mpr=squeeze(sum(reshape(T_mpr,mmax,[],mmax,2,2),2));

        %%
    %     T_mpl=(T_mpl+permute(T_mpl,[2 1 4 3]))/2;
%         T_mplr=(T_mpl)-mean(T_mpl(:));
% T_mplr=T_mpl;
    %     T_mprr=(T_mpr)-mean(T_mpr(:));

    end

    err=std(T_mplr(:))-std(T_ml(:));
    T_ml=T_mplr;
    T_mr=T_mprr;
    fprintf('error is %.4d ,corrf is %.4f, ratio is %.4f \n ',err,corrf,ra);    
    vsa=reshape(permute(T_ml,[1 3 2 4]),numel(T_ml)^0.5,[]);
    siz=size(vsa);
    vs=padarray(vsa,[2*mmax 2*mmax]-siz,0,'post');
    set(fi,'CData',real(vs));
    drawnow
    pause(0.5)
    end% T_ml=
    corrs=[corrs,corrf];
end
%%
imagesc(reshape(T_ml,numel(T_ml)^0.5,[]))
%%
imagesc(reshape(T_mpl,numel(T_ml)^0.5,[]))
