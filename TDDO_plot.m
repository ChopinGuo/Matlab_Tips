close all


figure(1);
plot(t,p(:,1),'r',t,p(:,2),'b',t,p(:,3),t,p(:,4),':',t,p(:,5),':','linewidth',2);
xlabel('time(s)');
ylabel('states and disturbance');
legend('x1','x2','x3','d1','d2');

% figure(2);
% plot(t,d(:,2),'r','linewidth',2);
