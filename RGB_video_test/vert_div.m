function G=vert_div(Z,nv)
%zero neumann boundary condition
%backward difference for vertical direction
G=reshape(Z,nv,[]);
G=reshape([zeros(1,numel(Z)/nv);G(2:end,:)-G(1:end-1,:)],size(Z));
end