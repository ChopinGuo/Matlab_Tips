tspan = [0,2];
y0 = [0,10];

[T,Y] = ode45('myfun_ode45',tspan,y0);

plot(T,Y(:,1)); % ��Y�ĵ�һ�����ݣ���y
hold on
plot(T,Y(:,2)); % ��Y�ĵ�һ�����ݣ���dy/dt