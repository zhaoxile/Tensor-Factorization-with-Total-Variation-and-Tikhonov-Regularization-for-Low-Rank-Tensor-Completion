function y=orthog_trans_2d(x,nv,nh)
y=reshape(dct(dct(reshape(x,[nv,nh,numel(x)/(nv*nh)]),[],1,'Type',3),[],2,'Type',3),size(x));
end