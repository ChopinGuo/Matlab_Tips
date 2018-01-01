% 关于ode45练习
% 网址：https://wenku.baidu.com/view/45a0a0b54b73f242326c5f7f.html

function Dx = myfun_ode45(t,x)

Dx = [0; 0];                  % 初始化Dx,使其为列数组
k = 1;

% 对微分方程的描述
Dx(1) = x(2);
Dx(2) = -k^2*x(1);