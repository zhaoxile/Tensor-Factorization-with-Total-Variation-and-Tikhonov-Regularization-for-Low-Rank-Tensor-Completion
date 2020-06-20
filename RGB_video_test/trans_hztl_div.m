function G=trans_hztl_div(Z,nv,nh)
%transpose operator of �hztl_div�
G=reshape(Z,[nv,nh,numel(Z)/(nv*nh)]);
G=reshape([-G(:,2,:),G(:,2:end-1,:)-G(:,3:end,:),G(:,end,:)],size(Z));
end