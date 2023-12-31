function [rho,i] = rho_k_thresholds(X,s1,s2,l,v_vec)

[~,N]=size(X);

[~,z]=size(v_vec);

k=sum(l);

n=N/k;

n_loop=1e3;

%The Hermit polynomial coefficients to compute p12 when abs(rho)<0.6
H=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 1;
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2, 0;
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,0, - 2;
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0, - 12,0;
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,16,0, - 48,0, 12;
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,32,0, - 160,0, 120,0;
0,0,0,0,0,0,0,0,0,0,0,0,0,0,64,0, - 480,0, 720,0, - 120;
0,0,0,0,0,0,0,0,0,0,0,0,0,128,0, - 1344,0, 3360,0, - 1680,0;
0,0,0,0,0,0,0,0,0,0,0,0,256,0, - 3584,0, 13440,0, - 13440,0, 1680;
0,0,0,0,0,0,0,0,0,0,0,512,0, - 9216,0, 48384,0, - 80640,0, 30240,0;
0,0,0,0,0,0,0,0,0,0,1024,0, - 23040,0,  161280,0, - 403200,0, 302400,0, - 30240;
0,0,0,0,0,0,0,0,0,2048,0, - 56320,0, 506880,0, - 1774080,0, 2217600,0, - 665280,0;
0,0,0,0,0,0,0,0,4096,0,- 135168,0, 1520640,0, - 7096320,0, 13305600,0, - 7983360,0, 665280;
0,0,0,0,0,0,0,8192,0, - 319488,0, 4392960,0, - 26357760,0, 69189120,0, - 69189120,0, 17297280,0;
0,0,0,0,0,0,16384,0, - 745472,0, 12300288,0, - 92252160,0, 322882560,0, - 484323840,0, 242161920,0, - 17297280;
0,0,0,0,0,32768,0, - 1720320,0, 33546240,0, - 307507200,0, 1383782400,0, - 2905943040,0, 2421619200,0, - 518918400,0;
0,0,0,0,65536,0, - 3932160,0, 89456640,0, - 984023040,0, 5535129600,0, - 15498362880,0, 19372953600,0, - 8302694400,0, 518918400;
0,0,0,131072,0,-8912896,0,233963520,0, -3041525760,0, 20910489600,0, -75277762560,0, 131736084480,0, - 94097203200,0, 17643225600,0;
0,0,262144,0,-20054016,0, 601620480,0,-9124577280,0, 75277762560,0, -338749931520,0, 790416506880,0, -846874828800,0, 317578060800,0, -17643225600;
0,524288,0, -44826624,0, 1524105216,0, -26671841280,0, 260050452480,0, - 1430277488640,0, 4290832465920,0, -6436248698880,0, 4022655436800,0, -670442572800,0;
1048576,0, - 99614720,0, 3810263040,0, -76205260800,0, 866834841600,0, - 5721109954560,0, 21454162329600,0, -42908324659200,0, 40226554368000,0, -13408851456000,0, 670442572800];

NH=20;

v1=v_vec/s1;
v2=v_vec/s2;

a1=v1/sqrt(2);
a2=v2/sqrt(2);

for i=1:z
    
a1_vec=fliplr(a1(i).^(0:NH))';
a2_vec=fliplr(a2(i).^(0:NH))';

a=(H*a1_vec).*(H*a2_vec)./2.^(1:NH+1)'./factorial(1:NH+1)';
a=4/pi*exp(-(a1(i)^2+a2(i)^2))*a;
     
A(:,i)=a;

end

Y=zeros(k,2,n);

for i=1:k
    Y(i,:,:)=X(:,(i-1)*n+1:i*n);
end

for i=1:k
    n12(i)=sum((Y(i,1,:)+1).*(Y(i,2,:)+1))/4;
    n1(i)=sum(Y(i,1,:)+1)/2;
    n2(i)=sum(Y(i,2,:)+1)/2;
end

for i=1:k
    c(i,:)=[n12(i),n1(i)-n12(i),n-n1(i)-n2(i)+n12(i),n2(i)-n12(i)];
end

T=zeros(z,k);
for i=1:z
    T(i,sum(l(1:i-1))+1:sum(l(1:i)))=1;
end

c=T*c;

p1=qfunc(v1);
p2=qfunc(v2);


u2=(2*p1-1).*(2*p2-1);


rho=0;

p12=zeros(1,z);


t=1;


for i=1:n_loop
    
    de=1/2/pi/sqrt(1-rho^2)*exp(-(v1.^2+v2.^2-2*rho*v1.*v2)/(1-rho^2)/2);
    
    if abs(rho)>0.6
        
        for j=1:z
            p12(j)=orthoscheme(1,1,rho,v1(j),v2(j));
        end
    else
        for j=1:z
            m12=rho.^(1:NH+1)*A(:,j)+u2(j);
            p12(j)=(m12-1+2*p1(j)+2*p2(j))/4;
        end
    end
    
    g=[de./p12;-de./(p1-p12);de./(1-p1-p2+p12);-de./(p2-p12)];
    
    de_sum=trace(c*g)/n;
    
    if abs(de_sum)<1e-4    %the precision can be adjusted to meet with demands
        
        break
    end
    
    
    de2=1/2/pi/sqrt(1-rho^2)*(rho/(1-rho^2)+(v1.*v2*(1-rho^2)-rho*(v1.^2-2*rho*v1.*v2+v2.^2))/(1-rho^2)^2).*exp(-(v1.^2+v2.^2-2*rho*v1.*v2)/(1-rho^2)/2);
    
    g2=[(de2.*p12-de.^2)./p12.^2;(-de2.*(p1-p12)-de.^2)./(p1-p12).^2;(de2.*(1-p1-p2+p12)-de.^2)./(1-p1-p2+p12).^2;(-de2.*(p2-p12)-de.^2)./(p2-p12).^2];
    
    
    de2_sum=trace(c*g2)/n;
    
    rho_0=rho-de_sum/de2_sum*t;
    
    while abs(rho_0)>1
        
        sr=0.5+0.1*rand;
        rho_0=sr*sign(rho_0)+(1-sr)*rho;
        
    end
    
    rho=rho_0;
    
end




end