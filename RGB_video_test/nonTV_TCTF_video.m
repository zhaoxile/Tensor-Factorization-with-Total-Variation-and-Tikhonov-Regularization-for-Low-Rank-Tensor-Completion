function [TC,psnrarr,ssimarr,iter]=nonTV_TCTF_video(X,Omega)
data=X(Omega);
known=Omega;
[nv,nh,nt,nb]=size(X);
ntvopts = [];
ntvopts.maxIter=500;
ntvopts.tol = 1e-5; 
ntvopts.alpha_adj = 0;
ntvopts.rank_adj = -ones(1,nb);
ntvopts.rank_inc = ones(1,nb);
ntvopts.rank_min = [27,ones(1,5),ones(1,max(nb-1-5,0))];
ntvopts.rank_max = 40*ones(1,nb);
Rankval=30;
[~,~,TC,~,~,iter] = old_TCTF_solver(data,known,[nv*nh,nt,nb],Rankval*ones(1,nb),ntvopts,X);
TC=reshape(TC,[nv,nh,nt,nb]);
psnrarr= TensorPSNR(TC,X);
ssimarr=TensorSSIM(TC,X);
end