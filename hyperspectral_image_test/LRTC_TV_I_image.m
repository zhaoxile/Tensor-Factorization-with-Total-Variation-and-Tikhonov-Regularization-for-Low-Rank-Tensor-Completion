function [TC,psnrtv1,ssimtv1]=LRTC_TV_I_image(Omega,V)
torder=3;
lambda=0.02;
lrtcalpha=[1/torder, 1/torder, 1/torder];
lrtcbeta=[1,1,0];
[TC]=LRTC_TV_I(Omega, V(Omega), lambda, lrtcalpha, lrtcbeta,size(V), torder );
psnrtv1=TensorPSNR(TC,V);
ssimtv1=TensorSSIM(TC,V);
end