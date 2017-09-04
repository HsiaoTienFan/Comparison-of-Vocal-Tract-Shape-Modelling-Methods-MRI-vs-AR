function spec = f_spec(s_r)

data_set = (s_r.^2).*pi;
data_set = data_set';
r_coef = calc_reflect_coef(data_set);
r_lpc = rc2lpc(r_coef);
spec = lpc2spec(r_lpc, 512);

