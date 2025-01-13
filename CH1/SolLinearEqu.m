% 07222113邓俊烽
% 解决线性方程组通用程序
function x = SolLinearEqu(A,b,method)
if nargin == 2
	method = 'ColLUGauss';
end
%nargin 是 MATLAB 中的一个内置变量，表示函数实际接收的输入参数的数量。这段代码检查 nargin 是否等于 2，即函数是否只接收了 2 个输入参数。
switch method
	case 'LUGauss'
		[L,U] = LUGauss(A)
		y = ForwardMethod(L,b);
		x = BackwardMethod(U,y);
	case 'FullLUGauss'
		[L,U,P,Q] = FullLUGauss(A);
		y = ForwardMethod(L,P*b);
		x = Q*BackwardMethod(U,y);
	case 'ColLUGauss'
		[L,U,P] = ColLUGauss(A)
		y = ForwardMethod(L,P*b);
		x = BackwardMethod(U,y);
	case 'Cholesky'
		L = Cholesky(A);
		y = ForwardMethod(L,b);
		x = BackwardMethod(L',y);
	case 'ImprovedCholesky'
		[L,D] = ImprovedCholesky(A);
		y = ForwardMethod(L,b);
		x = BackwardMethod(D*L',y);
	otherwise
		disp('算法不可用');
		disp('支持的为"LUGauss","FullLUGauss","ColLUGauss","Cholesky","ImprovedCholesky" ');
		x = 0;
end
