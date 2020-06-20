function [B,fftB]=B_ADMM_inn_solver(X,B,C,alpha,beta,rho,sigma,D_Sigma,time)
nv=size(B,1);
nr=size(B,2);
nt=size(B,3);
thresh=alpha/beta;
X=fft(X,[],3);
C=fft(C,[],3);
prerhs=zeros(nv,nr,nt);
Uc=zeros(nr,nr,nt);
Vc=zeros(nr,nt);
fftB=fft(B,[],3);
for ii=1:nt
   prerhs(:,:,ii)= X(:,:,ii)*(C(:,:,ii)')+rho*fftB(:,:,ii);
   [Uc(:,:,ii),temp]=svd((rho+sigma)*eye(nr)+C(:,:,ii)*(C(:,:,ii)'));
   Vc(:,ii)=diag(temp,0);
end
clear X C temp


Y=zeros(nv,nr,nt);
while time>0
    %% update P by soft thresholding
    P=vert_div(B,nv)-(1/beta)*Y;
    P=(abs(P)-thresh).*sign(P).*(abs(P)>thresh);
    
    
    %% update B
    % computing D_v^T(beta*P+Y) in Fourier domain
    RHS=fft(trans_vert_div(beta*P+Y,nv),[],3);
    
    for ii=1:1:nt
    %computing the right hand side
    F=prerhs(:,:,ii)+RHS(:,:,ii);
    %Transforming the right side by eigen-matrix of the linear system
    F=dct(F,[],1,'Type',2)*Uc(:,:,ii);
    %solving the diagonal linear system
    F=F./(repmat(D_Sigma,1,nr)+repmat(Vc(:,ii).',nv,1));
    %transform back to get the solution
    B(:,:,ii)=dct(F,[],1,'Type',3)*(Uc(:,:,ii)');
    end
    
    fftB=B;
    B=real(ifft(B,[],3));
    
    %% update Y and Z
    Y=Y+beta*(P-vert_div(B,nv));
    time=time-1;
end

end