function [C,fftC]=C_ADMM_inn_solver(X,B,C,sizeC,alpha,beta,rho,sigma,D_Sigma,time)
nr=sizeC(1);
nb=sizeC(2);
nt=sizeC(3);
thresh=alpha/beta;
X=fft(X,[],3);
B=fft(B,[],3);
prerhs=zeros(nr,nb,nt);
UB=zeros(nr,nr,nt);
VB=zeros(nr,nt);
fftC=fft(C,[],3);
for ii=1:nt
   prerhs(:,:,ii)= B(:,:,ii)'*X(:,:,ii)+rho*fftC(:,:,ii);
   [UB(:,:,ii),temp]=svd((rho+sigma)*eye(nr)+B(:,:,ii)'*B(:,:,ii));
   VB(:,ii)=diag(temp,0);
end
clear X B temp



Y=zeros(size(C));
while time>0
    %% update P
    P=spectr_div(C,nr,nb)-(1/beta)*Y;
    P=(abs(P)-thresh).*sign(P).*(abs(P)>thresh);
    
    
    %% update C
    % computing (beta*P+Y)D_s^T in Fourier domain
    RHS=fft(trans_spectr_div(beta*P+Y,nr,nb),[],3);
    
    for ii=1:1:nt
    %computing the right hand side
    F=prerhs(:,:,ii)+RHS(:,:,ii);
    %Transforming the right side by eigen-matrix of the linear system
    F=dct((UB(:,:,ii)'*F)',[],1,'Type',2)';
    %solving the diagonal linear system
    F=F./(repmat(VB(:,ii),1,nb)+repmat(D_Sigma,nr,1));
    %transform back to get the solution 
    C(:,:,ii)=dct((UB(:,:,ii)*F)',[],1,'Type',3)';
    end
    
    fftC=C;
    C=real(ifft(C,[],3));
    
    %% update Y
    Y=Y+beta*(P-spectr_div(C,nr,nb));
    time=time-1;
end

end