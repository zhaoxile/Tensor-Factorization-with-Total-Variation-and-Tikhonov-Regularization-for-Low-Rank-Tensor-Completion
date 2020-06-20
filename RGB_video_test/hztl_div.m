function G=hztl_div(Z,nv,nh)
%zero neumann boundary condition
%backward difference for horizontal direction
num=numel(Z);
G=reshape(Z,[nv,nh,num/(nv*nh)]);
sizg=size(G);
G=reshape([zeros(sizg(1),1,sizg(3)),G(:,2:end,:)-G(:,1:end-1,:)],size(Z));
end