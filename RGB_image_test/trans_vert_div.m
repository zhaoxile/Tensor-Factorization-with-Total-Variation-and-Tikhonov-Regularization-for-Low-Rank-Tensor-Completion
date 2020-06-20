function G=trans_vert_div(Z,~)
%transpose operator of vert_div
G=[-Z(2,:,:);Z(2:end-1,:,:)-Z(3:end,:,:);Z(end,:,:)];
end