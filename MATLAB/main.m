

a = importdata('M09FEXHA.000');


data_set_1x = a(1).data(: ,1);
data_set_1y = a(1).data(: ,2);

setappdata(0,'dataset_1x',data_set_1x)
setappdata(0,'dataset_1y',data_set_1y)

test = calc_reflect_coef(data_set_1y)