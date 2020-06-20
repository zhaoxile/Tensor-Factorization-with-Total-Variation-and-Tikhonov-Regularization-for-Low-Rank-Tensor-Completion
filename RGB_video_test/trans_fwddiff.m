function vout=trans_fwddiff(vin)
vout=[vin(1,:);vin(2:end,:)-vin(1:end-1,:);-vin(end,:)];
end