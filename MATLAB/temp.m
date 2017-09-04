Data = sel_file(datatype);


data_set_1x = ext_data(Data, datatype, 1);
data_set_1y = ext_data(Data, datatype, 2);



interp_1x = (data_set_1x(1)):0.01:(data_set_1x(size(data_set_1x,1)));
interp_1y = interp1(data_set_1x,data_set_1y,interp_1x,'spline');



length1 =  data_set_1x(size(data_set_1x,1));


tubelength1 = length1/11;
s1_r = get_tubes(interp_1x, interp_1y, tubelength1, start, finish, data_set_1x, datatype);




l_spec_1 = f_spec(s1_r);




function spec = f_spec(s_r)

data_set = s_r;
data_set = data_set';
r_coef = calc_reflect_coef(data_set);
r_coef = r_coef(1:size(r_coef,2)-1);
r_lpc = rc2lpc(r_coef);
spec = lpc2spec(r_lpc, 128);




function l_spec = lpc2spec(r_lpc, n)

order = 1:size(r_lpc,2);
tn = [0:(0.5/(n)):0.5];
zpow = 2i*pi*tn;
V1 = -order'*zpow;
V2 = V1(size(V1,1),:);
V3 = r_lpc'*exp(V2);
V4 = V3(1,:)*10;
hz = (1-exp(-zpow))./(1-V4);
l_spec = abs(hz);

