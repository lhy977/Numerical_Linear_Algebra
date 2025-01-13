% 07222113邓俊烽
% 解决最小二乘问题通用程序92
function x = SolLSQ(A,b)
[m,n] = size(A);
[Q,R] = QRDec(A);
disp('norm(Q*[R;zeros]-A)');
disp(norm(Q*[R;zeros(m-n,n)]-A));
if rank(R) < n
    disp('R 是奇异矩阵');
else
    disp('R 是非奇异矩阵');
end
Q1 = Q(:,1:n);
c1 = Q1'*b;
x = BackwardMethod(R,c1);
end
