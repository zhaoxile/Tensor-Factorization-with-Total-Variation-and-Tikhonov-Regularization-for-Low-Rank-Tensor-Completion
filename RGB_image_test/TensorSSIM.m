function ssimval=TensorSSIM(Xrecover,Xref)
[~,~,nb] = size(Xrecover);
ssimval=0;
for ii=1:nb
   ssimval=ssimval+ssim(Xrecover(:,:,ii),Xref(:,:,ii)); 
end
ssimval=ssimval/nb;
end