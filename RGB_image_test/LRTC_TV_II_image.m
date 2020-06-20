function [TC,psnrtv2,ssimtv2]=LRTC_TV_II_image(Omega,V)
N=3;
alpha=[1/N, 1/N, 1/N];
beta=[1,1,0];
lambda_1=0.5;
lambda_2=1000;
TC=LRTC_TV_II(Omega, V(Omega), lambda_1, lambda_2 ,alpha, beta, size(V), N );
psnrtv2=TensorPSNR(Z_TRLRTV2,V);
ssimtv2=TensorSSIM(Z_TRLRTV2,V);
end