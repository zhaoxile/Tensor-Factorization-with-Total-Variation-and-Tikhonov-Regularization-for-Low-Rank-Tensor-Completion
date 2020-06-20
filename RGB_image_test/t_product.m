function BC=t_product(fftB,fftC)
n1=size(fftB,1);
n2=size(fftC,2);
n3=size(fftB,3);
BC=zeros(n1,n2,n3);
for ii=1:n3
   BC(:,:,ii)=fftB(:,:,ii)*fftC(:,:,ii); 
end
BC=real(ifft(BC,[],3));
end