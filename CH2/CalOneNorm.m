%07222113邓俊烽
%算法2.5.1估计一范数
function Norm1 = CalOneNorm(A)

% A is a n*n matrix

[m,n] = size(A);
if m~=n
	disp('A不是n阶矩阵');
	Norm1 = 0;
	return;
end
times = 15; %迭代次数
temp = zeros(1,times);
for k = 1:times
	x = rand(n,1);
	x = x/sum(abs(x));%归一化处理使其一范数为1
	Target = 1;%内层循环的标志
	while Target == 1
		w = SolLinearEqu(A',x);%求解方程组A'w=x
		v = sign(w);
		z = SolLinearEqu(A,v);%求解方程组Az=v;
		if max(abs(z)) <= z'*x
			Target = 0;
			Norm1 = sum(abs(w));%退出迭代条件
		else%否则，更新 x 为z 中绝对值最大的分量对应的单位向量，并继续迭代。
			[~,maxindex] = max(abs(z));
			x = zeros(n,1);
			x(maxindex) = 1;
			Norm1 = sum(abs(w));
		end
	end
	temp(k) = Norm1;
end
Norm1 = max(temp);
end
