function r_coef = calc_reflect_coef(s_r)

reflect_coef(size(s_r,1)) = 0;
reflect_coef(1) = 1;

for i = 2:size(s_r,1)
    reflect_coef(i) = (s_r(i-1)-s_r(i))/(s_r(i-1)+s_r(i));
end
reflect_coef(size(s_r,1)+1) = 1;
r_coef = reflect_coef;

setappdata(0,'rcoef',r_coef);

