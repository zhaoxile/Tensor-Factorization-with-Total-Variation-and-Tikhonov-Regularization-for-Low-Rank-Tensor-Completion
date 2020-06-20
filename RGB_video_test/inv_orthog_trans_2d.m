function y=inv_orthog_trans_2d(x,nv,nh)
y=reshape(dct(dct(reshape(x,[nv,nh,numel(x)/(nv*nh)]),[],1),[],2),size(x));
end