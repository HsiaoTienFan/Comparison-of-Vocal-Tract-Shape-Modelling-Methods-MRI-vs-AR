function r_coef = calc_reflect_coef(data_set_1y)

reflect_coef(size(data_set_1y,1)) = 0;
reflect_coef(1) = 1;

for i = 2:size(data_set_1y,1)
    reflect_coef(i) = (data_set_1y(i)-data_set_1y(i-1))/(data_set_1y(i)+data_set_1y(i-1));
end

reflect_coef((size(reflect_coef,2))+1) = 1;
r_coef = reflect_coef;

setappdata(0,'rcoef',r_coef);

