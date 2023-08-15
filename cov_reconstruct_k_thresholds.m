function R = cov_reconstruct_k_thresholds(X,l,v_vec)

%Author: Yu-Hang Xiao
%Paper: Covariance Matrix Recovery from One-Bit Data with Non-Zero Quantization Thresholds: Algorithm and Performance Analysis
%Contact: yuhangxiao@szu.edu.cn

% l: the number of sub-intervals
% v_vec: the vector that contains the value of the threshold in each sub-interval
% X: the one-bit quantized data

[m,~]=size(X);

R=zeros(m,m);

s=zeros(1,m);

for i=1:length(l)
    v(sum(l(1:i-1))+1:sum(l(1:i)))=kron(v_vec(i),ones(1,l(i)));
end


for i=1:m
    
    s(i)=sigma_Newton_k_thresholds(X(i,:),v); % estimating the diagonal elementes
    
    R(i,i)=s(i)^2;
    
    for j=1:i-1
        rho=rho_k_thresholds([X(i,:);X(j,:)],s(i),s(j),l,v_vec);
        R(i,j)=rho*s(i)*s(j); % estimating the non-diagonal elementes
    end
end

R=R+R'-diag(diag(R));