%07222113邓俊烽
%利用隐式QR求解特征值，再用反幂法求解特征向量
function [eigvalue,eigvector] = FindEig(A)
n = length(A);

[Q,H]=HessenbergDec(A);
[~,~,kexi] = ImpQR(H);
eigvalue = zeros(n,1);
eigvector = zeros(n);
for i=1:n
    temp_H = H - eye(n)*kexi(i);%这里的最小特征值就是0了，就让它收敛到特征向量
    [eigvalue(i),eigvector(:,i)] = InversePowerMethod(temp_H);
end
eigvalue = eigvalue + kexi;
eigvector=Q*eigvector;
end