close all;

figure(1);
subplot(211)
plot(t,x(:,1),'r',t,x(:,4),'k:','linewidth',2);
% plot(t,x(:,1),'r','linewidth',2);
xlabel('time(s)');
ylabel('th and x1p');
legend('practical positiion signal', 'position signal estimation')
subplot(212)
plot(t,x(:,2),'r',t,x(:,5),'k:','linewidth',2);
xlabel('time(s)');
ylabel('dth and x2p');
legend('practical speed signal', 'speed sinal estimation')

figure(2);
plot(t,x(:,3),'k',t,x(:,6),'r','linewidth',2);
xlabel('time(s)');
ylabel('f and fp');
legend('practical uncertain part', 'uncertain part estimation')