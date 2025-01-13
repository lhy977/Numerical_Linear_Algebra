%07222113邓俊烽
%Page103GS迭代法
function x = GSIterate(A,b)
kmax = 1e7;%迭代次数
eps = 1e-18;%终止精度
tic%用于计时

n = length(A);
D = diag(diag(A));
L = - tril(A,-1);
U = - triu(A,1);
InvLU=inv(D-L);
A = InvLU*U;
a = InvLU*b;

k = 0;
x = zeros(n,1);
deltax = sqrt(eps^2*b'*b) + 1;
while deltax'*deltax > eps^2*(b'*b)
	k=k+1;
	if k>kmax
		disp(['已经达到最大迭代次数 ',num2str(kmax)]);
		break;
	end
	y = A*x + a;
	deltax = y - x;
	x = y;
end
disp(['迭代次数为',num2str(k)]);
disp(['运算时间为 ',num2str(toc),'s']);
end