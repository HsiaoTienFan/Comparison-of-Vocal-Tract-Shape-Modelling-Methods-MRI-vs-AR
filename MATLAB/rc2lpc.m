function r_lpc = rc2lpc(r_coef)

lpc_order = size(r_coef,2);

lin_pred_coef(1,1) = r_coef(1);


for i = 2:lpc_order
lin_pred_coef(i,i) = r_coef(i);
for j = 1:(i-1)
    lin_pred_coef(i,j) = lin_pred_coef(i-1,j)-(r_coef(i)*lin_pred_coef(i-1,i-j));
    
end
end

r_lpc = lin_pred_coef(lpc_order,:);


