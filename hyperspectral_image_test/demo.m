clear
addpath(genpath('test_msis'));
addpath(genpath('image_results'));
samplingrates=[0.05,0.1,0.15];
load('pompoms.mat');
imagnames=cell(1,1);
imagnames{1}='pompoms';
for kk=1:1
   
disp(imagnames{kk})
%% load data
V=V/max(abs(V(:)));
[nv,nh,nb]=size(V);
showchnlind=[1,2,31];
figure(1)
imshow(V(:,:,showchnlind),[])
hold on 
title('ground truth')
hold off

for jj=1:length(samplingrates)
p=samplingrates(jj);
disp('sampling rate:')
disp(p)
Omega = uint32(find(rand(numel(V),1)<p));
observe=zeros(size(V));
observe(Omega)=V(Omega);
figure(2);

imshow(observe(:,:,showchnlind),[])
hold on 
title(['Observation with p=',num2str(p)]);
hold off
clear observe

%% tensor completion
disp('Results of TCTF:')
disp('time cost of TCTF:')
tic
[TC,psnrsnontv,ssimsnontv]=nonTV_TCTF_video(V,Omega);
toc
figure(3);

imshow(TC(:,:,showchnlind),[])
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
[TC,psnrtnn,ssimtnn]=TNNmodel(V,Omega);
toc
figure(4);
imshow(TC(:,:,showchnlind),[])
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
[TC,psnrtvk,ssimtvk] = TV_TCTF_image(V,Omega);
toc
figure(5);
imshow(TC(:,:,showchnlind),[])
hold on 
title('Result of TCTF-TVK')
hold off
disp('psnr of TCTF-TVK:')
disp(psnrtvk)
disp('ssim of TCTF-TVK:')
disp(ssimtvk)



disp('LRTC-TV-I:')
disp('time cost of LRTC-TV-I:')
tic
[TC,psnrtv1,ssimtv1]=LRTC_TV_I_image(Omega,V);
toc
figure(6);
imshow(TC(:,:,showchnlind),[])
hold on
title('Result of LRTC-TV-I')
hold off

disp('psnr of LRTC-TV-I:')
disp(psnrtv1)
disp('ssim of LRTC-TV-I:')
disp(ssimtv1)

disp('LRTC-TV-II:')
disp('time cost of LRTC-TV-II:')
tic
[TC,psnrtv2,ssimtv2]=LRTC_TV_I_image(Omega,V);
toc
figure(7);
imshow(TC(:,:,showchnlind),[])
hold on
title('Result of LRTC-TV-II')
hold off

disp('psnr of LRTC-TV-II:')
disp(psnrtv2)
disp('ssim of LRTC-TV-II:')
disp(ssimtv2)













end
end


