%07222113邓俊烽
%Page185带原点位移的QR
function eigvalues = OQRIterate(A)
kmax = 1e7;
eps = 1e-12;
tic; 

n = size(A, 1);
A_k = A; 
for k = 1:kmax
    % 原点位移：取 A_k 的右下角元素作为位移量
    mu = A_k(end, end); % 位移量
    [Q, R] = QRDec(A_k - mu * eye(n)); 
    A_k = R * Q + mu * eye(n);
    
    lower = tril(A_k, -1); 
    if max(abs(lower(:))) < eps
        break; 
    end
end
eigvalues = diag(A_k);
disp(['迭代次数为 ', num2str(k)]);
disp(['迭代时间为 ', num2str(toc), 's']);
end