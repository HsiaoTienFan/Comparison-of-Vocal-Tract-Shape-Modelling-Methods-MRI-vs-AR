function test = rc2lpc(r_coef)
x = r_coef'
x = x(1:size(x,1)-1)
hrc2lpc = dsp.RCToLPC;
[A] = step(hrc2lpc, x);
A = A(2:size(A))
test = A';


