function G=trans_hztl_div(Z,~,~)
%transpose of "hztl_div" operator
G=[-Z(:,2,:),Z(:,2:end-1,:)-Z(:,3:end,:),Z(:,end,:)];
end