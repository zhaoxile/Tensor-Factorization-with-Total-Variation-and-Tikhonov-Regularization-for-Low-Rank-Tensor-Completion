function [X,psnrarr,ssimarr,iter]=TV_TCTF_solver(data,known,Dims,Rankval,opts)

nv=Dims(1);
nh=Dims(2);
nb=Dims(3);
btime=1;
ctime=1;
thetav=pi/nv;
sigmav=(2*opts.beta)*(1-cos(((1:nv).'-1)*thetav));
thetah=pi/nh;
sigmah=(2*opts.beta)*(1-cos(((1:nh)-1)*thetah));
psnrarr=zeros(1,opts.maxIter);
ssimarr=zeros(1,opts.maxIter);
% initialization of B, C and X^k
% ntvopts = [];
% ntvopts.maxIter=20;
% ntvopts.tol = -1e-5; 
% ntvopts.alpha_adj = 0;
% ntvopts.rank_adj = -ones(1,nt);
% ntvopts.rank_inc = ones(1,nt);
% ntvopts.rank_min = [27,ones(1,5),ones(1,nt-1-5)];
% ntvopts.rank_max = 40*ones(1,nt);
% [temp1,temp2,X] = old_TCTF_solver(data,known,[nv*nh,nb,nt],Rankval*ones(1,nt),ntvopts);
% X(known)=data;
% B=zeros(nv*nh,Rankval,nt);
% C=zeros(Rankval,nb,nt);
% for ii=1:nt
%   B(:,:,ii)=temp1{ii};
%   C(:,:,ii)=temp2{ii};
% end
%  
% clear temp1 temp2
X=rand(nv,nh,nb);
X(known)=data;
temp=fft(X,[],3);
B=zeros(nv,Rankval,nb);
C=zeros(Rankval,nh,nb);
for ii=1:nb
   [U,S,V]=svds(temp(:,:,ii),Rankval);
   S=diag(S,0);
   B(:,:,ii)=U.*repmat(S.',nv,1);
   C(:,:,ii)=V';
end
B=real(ifft(B,[],3));
C=real(ifft(C,[],3));
clear temp U S V
   % iter=0;
% 
% %% The outer iteration
%while 1
  for ii=1:opts.maxIter  
   % B subproblem
   temp=B;
   [B,fftB]=B_ADMM_inn_solver(X,B,C,opts.alpha_phys,opts.beta,opts.rho,opts.sigma,sigmav,btime);
   curres=norm(B(:)-temp(:),2).^2/(norm(B(:),2).^2);
   % C subproblem
   temp=C;
   [C,fftC]=C_ADMM_inn_solver(X,B,C,opts.alpha_phys,opts.beta,opts.rho,opts.sigma,sigmah,ctime);
   curres=max(curres,norm(C(:)-temp(:),2).^2/(norm(C(:),2)).^2);
   % X subproblem
   temp=X;
   X=max((t_product(fftB,fftC)+opts.rho*X)/(1+opts.rho),0);
     X(X>1)=1;
   X(known)=data;
%    psnrarr(ii)= TensorPSNR(X,Xref);
%    ssimarr(ii)=TensorSSIM(X,Xref);
   curres=max(curres,norm(X(:)-temp(:),2).^2/(norm(X(:),2)).^2);
    %iter=iter+1;

% B subproblem
   
%    [B,fftB]=B_ADMM_inn_solver(X,B,C,[nv,nh,Rankval,nt],opts.alpha_phys,opts.beta,opts.rho,Sigma_2d,btime);
%   
%    % C subproblem
%    
%    [C,fftC]=C_ADMM_inn_solver(X,B,C,[Rankval,nb,nt],opts.alpha_spec,opts.beta,opts.rho,Sigma_s,ctime);
%    
%    % X subproblem
%    %temp=X;
%    
%    X=max((t_product(fftB,fftC)+opts.rho*X)/(1+opts.rho),0);
%    %X=max((t_product(fftB,fftC)+opts.rho*X),0);
%    X=X/max(abs(X(:)));
% curres=norm(X(known)-data)/normdata
   if curres<=opts.tol
      break; 
   end
  end
 iter=ii;
end