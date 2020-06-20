function vout=fwddiff(vin)
vout=vin(1:end-1,:)-vin(2:end,:);
end