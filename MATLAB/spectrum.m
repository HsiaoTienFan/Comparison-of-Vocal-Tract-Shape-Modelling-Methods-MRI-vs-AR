function r_coef = calc_reflect_coef(data_set_1y)

reflect_coef(size(data_set_1y,1)) = 0;
reflect_coef(1) = 1;

for i = 2:size(data_set_1y,1)
    reflect_coef(i) = (data_set_1y(i-1)-data_set_1y(i))/(data_set_1y(i-1)+data_set_1y(i));
end

reflect_coef((size(reflect_coef,2))+1) = 1;
r_coef = reflect_coef;

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

function l_spec = lpc2spec(r_lpc, n)

order = 1:size(r_lpc,2);
tn = [0:(0.5/(n-1)):0.5];
zpow = 2i*pi*tn;
hz = (1-exp(-zpow))./(1-(r_lpc*exp(-order'*zpow)));

l_spec = abs(hz);


function setup()

fs = 34000*11/(2*17.59);

x = (1:n)*fs/(2*n);


plot (x, log(l_spec))

