% ����ode45��ϰ
% ��ַ��https://wenku.baidu.com/view/45a0a0b54b73f242326c5f7f.html

function Dx = myfun_ode45(t,x)

Dx = [0; 0];                  % ��ʼ��Dx,ʹ��Ϊ������
k = 1;

% ��΢�ַ��̵�����
Dx(1) = x(2);
Dx(2) = -k^2*x(1);