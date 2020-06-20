function [B,fftB]=B_ADMM_inn_solver(X,B,C,sizeB,alpha,beta,rho,sigma,D_Sigma,time)
nv=sizeB(1);
nh=sizeB(2);
nr=sizeB(3);
nt=sizeB(4);
thresh=alpha/beta;
X=fft(X,[],3);
C=fft(C,[],3);
prerhs=zeros(nv*nh,nr,nt);
Uc=zeros(nr,nr,nt);
Vc=zeros(nr,nt);
fftB=fft(B,[],3);
for ii=1:nt
   prerhs(:,:,ii)= X(:,:,ii)*(C(:,:,ii)')+rho*fftB(:,:,ii);
   [Uc(:,:,ii),temp]=svd((rho+sigma)*eye(nr)+C(:,:,ii)*(C(:,:,ii)'));
   Vc(:,ii)=diag(temp,0);
   %[Uc(:,:,ii),Vc(:,ii)]=svd(rho*eye(nr)+C(:,:,ii)*(C(:,:,ii)'));
end
clear X C temp



Y=zeros(size(B));
Z=zeros(size(B));
while time>0
    %% update P and Q by soft thresholding
    P=vert_div(B,nv)-(1/beta)*Y;
    P=(abs(P)-thresh).*sign(P).*(abs(P)>thresh);
    Q=hztl_div(B,nv,nh)-(1/beta)*Z;
    Q=(abs(Q)-thresh).*sign(Q).*(abs(Q)>thresh);
%     Q=vert_div(B,nv)+hztl_div(B,nv,nh)-(1/beta)*Y;
%     Q=(abs(Q)-thresh).*sign(Q).*(abs(Q)>thresh);
    
    %% update B
    % computing D_1^T(beta*P+Y) and D_2^T(beta*Q+Z) in Fourier domain
    RHS=fft(trans_vert_div(beta*P+Y,nv)+trans_hztl_div(beta*Q+Z,nv,nh),[],3);
    %RHS=fft(trans_vert_div(beta*Q+Y,nv)+trans_hztl_div(beta*Q+Y,nv,nh),[],3);
    for ii=1:1:nt
    %computing the right hand side
    F=prerhs(:,:,ii)+RHS(:,:,ii);
    %Transforming the right side by eigen-matrix of the linear system  
    F=reshape(dct(dct(reshape(F,[nv,nh,nr]),[],1),[],2),[nv*nh,nr])*Uc(:,:,ii);
    %solving the diagonal linear system
    F=F./(repmat(D_Sigma,1,nr)+repmat(Vc(:,ii).',nv*nh,1));
    %transform back to get the solution
    B(:,:,ii)=reshape(dct(dct(reshape(F,[nv,nh,nr]),[],1,'Type',3),[],2,'Type',3),[nv*nh,nr])*(Uc(:,:,ii)');
    end
    
    fftB=B;
    B=real(ifft(B,[],3));
    
    %% update Y and Z
    Y=Y+beta*(P-vert_div(B,nv));
    Z=Z+beta*(Q-hztl_div(B,nv,nh));
    time=time-1;
end

end