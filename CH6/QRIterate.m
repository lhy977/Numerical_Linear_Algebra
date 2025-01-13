%07222113邓俊烽
%Page174QR迭代
function eigvalues = QRIterate(A)
kmax = 1e7;
eps = 1e-12;
tic

n = size(A, 1);
A_k = A; 

for kmax = 1:kmax
    [Q, R] = QRDec(A_k); 
    A_k = R * Q;
    lower = tril(A_k, -1); 
    if max(abs(lower(:))) < eps
        break;
    end
end
eigvalues = diag(A_k);
disp(['迭代次数为',num2str(k)]);
disp(['迭代时间为',num2str(toc),'s']);
end