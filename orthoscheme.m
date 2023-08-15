function p = orthoscheme(s1,s2,s12,v1,v2)
% using numerical integration to compute the orthoscheme probability
w1=s1^2*s2^2-s12^2;

syms x1; syms x2

fun = @(x1,x2) 1/2/pi/sqrt(w1)*exp(-1/2/w1*(s2^2*x1.^2-2*s12*x1.*x2+s1^2*x2.^2));            
               
p=integral2(fun,v1,10,v2,10);    