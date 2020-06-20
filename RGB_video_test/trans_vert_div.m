function G=trans_vert_div(Z,nv)
%transpose operator of �vert_div�
G=reshape(Z,nv,[]);
G=reshape([-G(2,:);G(2:end-1,:)-G(3:end,:);G(end,:)],size(Z));
end