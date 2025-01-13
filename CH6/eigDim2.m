% 07222113邓俊烽
% 求一个二维矩阵的特征值
function [l1, l2] = eigDim2(A)
a = A(1, 1);
b = A(1, 2);
c = A(2, 1);
d = A(2, 2);
delta = sqrt((a + d)^2 - 4*(a*d - b*c));
l1 = (a + d + delta)/2;
l2 = (a + d - delta)/2;
end