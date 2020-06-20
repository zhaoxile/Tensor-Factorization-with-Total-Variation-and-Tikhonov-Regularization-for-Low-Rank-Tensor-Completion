function [ Xapp,psnrarr,ssimarr,iter] =TV_TCTF_image(X,Omega)

%% produce data

data=X(Omega);
known=Omega;
[nv,nh,nt,nb]=size(X);

%% our method 
opts = [];
opts.maxIter=500;
opts.tol = 1e-6; 

opts.alpha_phys = 5;
opts.alpha_spec = 5;
opts.rho= 5e-6;
opts.beta=1.5;
opts.sigma=0.04;


tRank=30;
Dims=[nv,nh,nt,nb];

[Xapp,~,~,iter] = TV_TCTF_solver(data,known,Dims,tRank,opts);
%% compute PSNR
Xapp=reshape(Xapp,Dims);
 psnrarr= TensorPSNR(Xapp,X);
 ssimarr=TensorSSIM(Xapp,X);
end



