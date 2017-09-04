
function setup(n, l_spec, x)

fs = 34000*11/(2*17.59);

x = (1:n)*fs/(2*n);


plot (x, log(l_spec))
