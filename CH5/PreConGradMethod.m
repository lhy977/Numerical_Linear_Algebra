%07222113邓俊烽
%算法5.4.1.预优共轭梯度法
function x = PreConGradMethod(A,b)
%A是一个对称矩阵
kmax = 1e7;
eps = 1e-40;
tic

n = length(b);
x = zeros(n,1); % initial value x0
r = b - A*x;
k = 0;

while r'*r > eps^2
	k=k+1;
	if k>kmax
		disp(['超出迭代次数',num2str(kmax)]);
		break;
    end
	z = r./diag(A);

	if k==1
		p = z;
		rz = r'*z;
	else
		temprz = rz;
		rz = r'*z;
		beta = rz/temprz;
		p = z + beta*p;
	end
	Ap = A*p;
	alpha = rz/(p'*Ap);
	x = x + alpha*p;
	r = r - alpha*Ap;
end
disp(['迭代次数为',num2str(k)]);
disp(['迭代时间为 ',num2str(toc),'s']);
end