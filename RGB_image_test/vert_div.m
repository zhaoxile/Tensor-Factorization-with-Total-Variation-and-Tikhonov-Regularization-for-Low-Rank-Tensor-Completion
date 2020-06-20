function G=vert_div(Z,~)
%zero neumann boundary condition
%backward difference for vertical direction
sizez=size(Z);
G=[zeros(1,sizez(2),sizez(3));Z(2:end,:,:)-Z(1:end-1,:,:)];
end