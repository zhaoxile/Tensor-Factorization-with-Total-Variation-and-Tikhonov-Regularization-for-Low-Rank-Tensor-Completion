clear
addpath(genpath('test_video'));
addpath(genpath('RGBvideoresults'));
samplingrates=[0.05,0.1,0.15];
videoname='salesman_qcif';
%% load data
load([videoname,'.mat']);
V=double(X(:,:,:,1:40));
clear X;
  
disp(videoname)

V=V/max(abs(V(:)));
[nv,nh,nb,nt]=size(V);
shownframe=ceil(nt/2);

figure(1)
imshow(V(:,:,:,shownframe),[])
hold on 
title('ground truth')
hold off
Vper=permute(V,[1,2,4,3]);
for jj=1:length(samplingrates)
p=samplingrates(jj);
disp('sampling rate:')
disp(p)
Omega = uint32(find(rand(numel(Vper),1)<p));
observe=zeros(size(Vper));
observe(Omega)=Vper(Omega);
observe=permute(observe,[1,2,4,3]);
figure(2);
imshow(observe(:,:,:,shownframe),[])
hold on 
title(['Observation with p=',num2str(p)]);
hold off
clear observe

%% tensor completion
disp('Results of TCTF:')
disp('time cost of TCTF:')
tic
[TC,psnrsnontv,ssimsnontv]=nonTV_TCTF_video(Vper,Omega);
toc
TC=permute(reshape(TC,[nv,nh,nt,nb]),[1,2,4,3]);
figure(3);

imshow(TC(:,:,:,shownframe),[])
hold on 
title('Result of TCTF')
hold off
disp('psnr of TCTF:')
disp(psnrsnontv)
disp('ssim of TCTF:')
disp(ssimsnontv)

disp('Result of TNN:')
disp('time cost of TNN:')
tic
[TC,psnrtnn,ssimtnn]=TNNmodel(Vper,Omega);
toc
TC=permute(reshape(TC,[nv,nh,nt,nb]),[1,2,4,3]);
figure(4);
imshow(TC(:,:,:,shownframe),[])
hold on 
title('Result of TNN')
hold off
disp('psnr of TNN:')
disp(psnrtnn)
disp('ssim of TNN:')
disp(ssimtnn)


disp('TCTF-TVK:')
disp('time cost of TCTF-TVK:')
tic
[TC,psnrtvk,ssimtvk] = TV_TCTF_image(Vper,Omega);
toc
TC=permute(reshape(TC,[nv,nh,nt,nb]),[1,2,4,3]);
figure(5);
imshow(TC(:,:,:,shownframe),[])
hold on 
title('Result of TCTF-TVK')
hold off
disp('psnr of TCTF-TVK:')
disp(psnrtvk)
disp('ssim of TCTF-TVK:')
disp(ssimtvk)



disp('LRTC-TV-I:')
disp('Caution: LRTC-TV-I is time-consuming for RGB video completion:')
disp('time cost of LRTC-TV-I:')
tic
[TC,psnrtv1,ssimtv1]=LRTC_TV_I_video(Omega,Vper);
toc
TC=permute(reshape(TC,[nv,nh,nt,nb]),[1,2,4,3]);
figure(6);
imshow(TC(:,:,:,shownframe),[])
hold on
title('Result of LRTC-TV-I')
hold off

disp('psnr of LRTC-TV-I:')
disp(psnrtv1)
disp('ssim of LRTC-TV-I:')
disp(ssimtv1)

disp('LRTC-TV-II:')
disp('Caution: LRTC-TV-II is very time-consuming for RGB video completion:')
disp('time cost of LRTC-TV-II:')
tic
[TC,psnrtv2,ssimtv2]=LRTC_TV_II_video(Omega,Vper);
toc
TC=permute(reshape(TC,[nv,nh,nt,nb]),[1,2,4,3]);
figure(7);
imshow(TC(:,:,:,shownframe),[])
hold on
title('Result of LRTC-TV-II')
hold off

disp('psnr of LRTC-TV-II:')
disp(psnrtv2)
disp('ssim of LRTC-TV-II:')
disp(ssimtv2)













end



