function [TC,psnrtnn,ssimtnn]=TNNmodel(V,Omega)
opts=[];
opts.tol = 1e-5;
opts.maxit =500;
opts.rho = 1.1;
opts.beta = 1e-2;
opts.max_beta = 1e10;
opts.Xtrue=V;
TC=rand(size(V));
TC(Omega)=V(Omega);
[TC,Out]=LRTC_tnn(TC,Omega,opts);
psnrtnn=Out.PSNR;
ssimtnn=Out.SSIM;
end