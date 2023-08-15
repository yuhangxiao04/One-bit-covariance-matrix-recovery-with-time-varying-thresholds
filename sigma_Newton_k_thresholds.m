function [s,i] = sigma_Newton_k_thresholds(x,v)

n_loop=1e3;

[~,k]=size(v);

[~,N]=size(x);

n=N/k;

Y=zeros(k,n);

for i=1:k
    Y(i,:)=x((i-1)*n+1:i*n);
end


for i=1:k
    
    np(i)=sum(Y(i,:)+1)/2;
    
end

p=np/n;

c=[p',1-p'];

o=1;s=0;
for i=1:k
    if p(i)<0.5
        s=s+v(i)/qfuncinv(p(i));
        o=o+1;
    end
end
s=s/o;

t=1;

for i=1:n_loop
    
    s0=s;
    
    p1=qfunc(v./s)';
    
    de=1/sqrt(2*pi)*v/s^2.*exp(-v.^2/2/s^2);
    
    g=[de'./p1,-de'./(1-p1)];
    
    de_sum=trace(g'*c);
    
    if abs(de_sum)<1e-4   %the precision can be adjusted to meet with demands
        
        break
    end
    
    
    de2=1/sqrt(2*pi)*exp(-v.^2/2/s^2).*(v.^3/s^5-2*v/s^3);
    
    
    g2=[(de2'.*p1-de'.^2)./p1.^2,(-de2'.*(1-p1)-de'.^2)./(1-p1).^2];
    
    de2_sum=trace(g2*c');
    
    s=s-de_sum/de2_sum*t;
    
    while s<0
        s=(s+s0)/2;
    end
    
    
end
