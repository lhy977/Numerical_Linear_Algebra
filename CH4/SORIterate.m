%07222113邓俊烽
%Page122SOR迭代法
function x = SORIterate(A,b,omega)
kmax = 1e5;
eps = 1e-18;
tic
n = length(A);
D = diag(diag(A));
L = - tril(A,-1);
U = - triu(A,1);
InverseDOL=inv(D-omega*L);
A = InverseDOL*((1-omega)*D + omega*U);
a = omega*InverseDOL*b;

k = 0;
x = zeros(n,1);
deltax = sqrt(eps^2*b'*b) + 1;
while deltax'*deltax > eps^2*(b'*b)
	k=k+1; % The number of iterations
	if k>kmax
		disp(['迭代达到最大次数 ',num2str(kmax)]);
		break;
	end
	y = A*x + a;
	deltax = y - x;
	x = y;
end
disp(['迭代次数为',num2str(k)]);
disp(['迭代时间为',num2str(toc),'s']);
end