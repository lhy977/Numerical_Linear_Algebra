%07222113邓俊烽
%算法5.1.1.最速下降法
function x = DesMethod(A,b)
%A是一个对称矩阵
kmax = 1e7;
eps = 1e-12;
tic

n = length(b);
x = zeros(n,1); 
r = b - A*x;
k = 0;
rr = r'*r;
while rr > eps^2*b'*b
	k=k+1;
	if k>kmax
		disp(['超出迭代次数',num2str(kmax)]);
		break;
	end
	Ar = A*r;
	alpha = rr/(r'*Ar);
	x = x + alpha*r;
	r = r - alpha*Ar;
	rr = r'*r;
end

disp(['迭代次数为',num2str(k)]);
disp(['迭代时间为 ',num2str(toc),'s']);
end