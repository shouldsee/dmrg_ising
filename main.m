K=0.44;
ss=[-1 1]
mmax=30;
global oper opts m
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
oper=@times;
T_ml=W;
figure(2)
fi=imagesc(zeros(2*mmax,2*mmax));
corrs=[];
% for i=1:1024;
% Ts=0.5:0.1:3;
Ks=-0.6:0.05:0.6;
% Ts=1./Ks;
% Ks=-0.20;
for K=Ks;
T=1/K;
err=1;

logW=-0.5*K*(msh1.*msh2+msh1p.*msh2p+msh1.*msh1p+msh2.*msh2p);
W=exp(logW);
T_ml=W;

while abs(err)>1e-15 | 2*m<mmax;
    m=size(T_ml,1);
    if 4*m<=2*mmax
        % tp=bsxfun(@times,T_ml,shiftdim(W,-2));
        tp=bsxfun(oper,T_ml,shiftdim(W,-2));
        tp=permute(tp,[1 3 2 4 5 6]);
        T_mplr=reshape(tp,2*m,2*m,2,2);
        T_mprr=T_mplr;
        
        [phi,psi]=superblock(T_ml,W);
        rhoL=squeeze(sum(sum(oper(phi, psi),3),4));
        corrf=NNcorr(phi,psi,ss);
    
    elseif opt.renorm==1;
        %% calculate T_2m
        %     tp=bsxfun(oper,bsxfun(oper,T_ml,shiftdim(W,-2)),permute(T_ml,[5:8 flip(1:4)]));
        %     tp=T_ml+shiftdim(W,-2)+permute(T_ml,[5:8 3 4 1 2]);
        %     figure(3)
            [phi,psi]=superblock(T_ml,W);
            rhoL=squeeze(sum(sum(oper(phi, psi),3),4));       
            corrf=NNcorr(phi,psi,ss);

            D=rot90(D,2);

        if strcmp(opt.decomp,'svd');

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
        T_mpl=oper(...
                    oper(...
                    oper(Ol,...
                        T_mlo),...
                        Ql ),...
                        permute((W),[5,6,1,7,2,8,3,4])...
                        );
        T_mpl=squeeze(sum(reshape(T_mpl,mmax,[],mmax,2,2),2));
        T_mpl=0.5*(T_mpl+permute(T_mpl,[2 1 4 3]));
        T_mplr=T_mpl;
    %     T_mplr=T_mpl-mean(T_mpl(:));

    else 
        err=0;
    end

    err=std(T_mplr(:))-std(T_ml(:));
    T_ml=T_mplr;
    T_mr=T_mprr;
    fprintf('error is %.4d ,corrf is %.4f, ratio is %.4f \n ',err,corrf,ra);    
    vsa=reshape(permute(T_ml,[1 3 2 4]),numel(T_ml)^0.5,[]);
    siz=size(vsa);
    vs=padarray(vsa,[2*mmax 2*mmax]-siz,0,'post');
%     set(fi,'CData',real(vs));
%     drawnow
%     pause(0.5)
    end% T_ml=
    corrs=[corrs,corrf];
end
%%
% plot(Ks,[0 diff(2*corrs)]);
plot(Ks,corrs)
% imagesc(reshape(T_ml,numel(T_ml)^0.5,[]))
%%
% imagesc(reshape(T_mpl,numel(T_ml)^0.5,[]))
