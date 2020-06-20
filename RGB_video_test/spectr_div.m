function G=spectr_div(Z,~,~)
%zero neumann boundary condition
%backward difference for spectral direction

sizez=size(Z);
G=[zeros(sizez(1),1,sizez(3)),Z(:,2:end,:)-Z(:,1:end-1,:)];
end