tspan = [0,2];
y0 = [0,10];

[T,Y] = ode45('myfun_ode45',tspan,y0);

plot(T,Y(:,1)); % 画Y的第一列数据，即y
hold on
plot(T,Y(:,2)); % 画Y的第一列数据，即dy/dt