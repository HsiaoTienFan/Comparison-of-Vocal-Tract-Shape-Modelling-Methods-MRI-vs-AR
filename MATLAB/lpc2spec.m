function l_spec = lpc2spec(r_lpc, n)

order = 1:size(r_lpc,2);
tn = [0:(0.5/(n)):0.5];
tn = tn(1:(size(tn,2)-1));
zpow = 2i*pi*tn;
hz = (1-exp(-zpow))./(1-(r_lpc*exp(-order'*zpow)));
l_spec = abs(hz);
