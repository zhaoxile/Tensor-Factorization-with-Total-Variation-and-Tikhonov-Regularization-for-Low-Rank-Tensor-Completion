function psnr = TensorPSNR(Xrecover,Xfull)

% Xrecover = max(0,Xrecover);
% Xrecover = min(maxP,Xrecover);
[nv,nh,nb] = size(Xrecover);
% MSE = norm(Xfull(:)-Xrecover(:))^2/(3*m*n*dim);
% psnr = 10*log10(maxP^2/MSE);
nim=nv*nh;
Xrecover=reshape(Xrecover,nim,[]);
Xfull=reshape(Xfull,nim,[]);
MSE=sum((Xrecover-Xfull).^2,1);
Maxpxs=nim*(max(abs(Xfull)).^2);
psnr=sum(10*log10(Maxpxs./MSE))/nb;
end