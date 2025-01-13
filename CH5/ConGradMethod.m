%07222113邓俊烽
%算法5.3.1.实用的共轭梯度法
function x = ConGradMethod(A,b)%conjugate gradient method
%A是一个对称矩阵
kmax = 1e7;
eps = 1e-30;
tic

n = length(b);
x = zeros(n,1);
r = b - A*x;
rr = r'*r;
k = 0;

while rr > eps^2*(b'*b)
	k=k+1;
	if k>kmax
		disp(['超出迭代次数',num2str(kmax)]);
		break;
	end
	if k==1
		p = r;
	else
		beta = rr/temprr;
		p = r + beta*p;
	end
	Ap = A*p;
	alpha = rr/(p'*Ap);
	x = x + alpha*p;
	r = r - alpha*Ap;
	temprr = rr;
	rr = r'*r;
end

disp(['迭代次数为',num2str(k)]);
disp(['迭代时间为',num2str(toc),'s']);
end
