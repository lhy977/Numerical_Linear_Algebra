%07222113邓俊烽
%Page166幂法
function [kxi,u] = PowerMethod(A)
%求解最大的特征值
kmax = 1e7;
eps = 1e-12;
tic

n = length(A);
u = ones(n,1);
k = 0;
while true
    k = k+1;
    if k>kmax
		disp(['超出迭代次数 ',num2str(kmax)]);
		break;
    end
    tempu=u;
    y=A*u;
    [~,pos]=max(abs(y));
    kxi=y(pos);
    u=y/kxi;
    if norm(u-tempu)<eps
        break;
    end
end
disp(['迭代次数为',num2str(k)]);
disp(['迭代时间为',num2str(toc),'s']);
end
