function [ ] = test_performance()
    filename=cell(8,1);
    filename{1}='TestImages/cutFlower.jpg';
    filename{2}='Burano.jpg';
    filename{3}='barbara';
    filename{4}='facade';
    filename{5}='house';
    filename{6}='lena';
    filename{7}='peppers';
    filename{8}='sailboat';
    
    for ff=1:1
        myName=filename{ff};
        disp(myName)
        A=imread(myName);
        figure(1);  imshow(A);

        myrate=0.95:-0.05:0.85;
        myResult=cell(2,numel(myrate));
        A=double(A)/255.0;

        for iterate=1:numel(myrate)
            rate=1 - myrate(iterate);
            [row, col, channel]=size(A);

            
            index = uint32(find(rand(numel(A),1)<rate));
            value=A(index);

            tsize=[row, col, channel];
            N=3;
            lambda=0.02;
            alpha=[1/N, 1/N, 1/N];
            beta=[1,1,0];

            fprintf('------TR_LRTV---------- \n');
            Z_TRLRTV=LRTC_TV_I(index, value, lambda, alpha, beta, tsize, N );
            disp('LRTC_TVI:')
            psnr=TensorPSNR(Z_TRLRTV,A)
            SSIM=TensorSSIM(Z_TRLRTV,A)
            figure(3);imshow(Z_TRLRTV);

            fprintf('--------------TR_LRTV2-------------------\n');
            lambda_1=0.5;
            lambda_2=1000;

            Z_TRLRTV2=LRTC_TV_II(index, value, lambda_1, lambda_2 ,alpha, beta, tsize, N );
            disp('LRTC_TVII:')
            psnr=TensorPSNR(Z_TRLRTV2,A)
            SSIM=TensorSSIM(Z_TRLRTV2,A)

%             myResult{1,iterate}=Z_TRLRTV;  %LRTC_TV_I
%             myResult{2,iterate}=Z_TRLRTV2; %LRTC_TV_II
% 
%             myName=['Result/' filename{ff} '_Result.mat'];
%             save(myName,'myResult','myrate');
% 
%             myName=['Result/' filename{ff}, '_Data.mat'];
%             save(myName,'A');
        end
    end
    
%     [perform_RSE,perform_PSNR] = performance_eval();
%     figure(1);plot(perform_RSE);
%     figure(2);plot(perform_PSNR);
end

