close all;
clear;
clc;



w0 = 1;
kappa = 1;

alpha1 = 3;
alpha2 = 3;
alpha3 = 1;

A0 = [0 1;0 0];
B0 = [0; 1];
G0_star = [w0*alpha1; w0^2*alpha2];
C0 = [1, 0];

[Ord Ord] = size(A0);


%% solve ARE

setlmis([]);
P = lmivar(1,[2 1]);
lmiterm([1 1 1 P],A0,1,'s');                    % LMI #1: A0*P+P*A0'
lmiterm([1 1 1 0],-2*(B0*B0'));                 % LMI #1: B0'  
lmiterm([-2,1,1,P],1,1);
lmis=getlmis;

[tmin, feas] = feasp(lmis);
PP=dec2mat(lmis,feas,P);
PP_inv = inv(PP);
K = -B0'*PP_inv;

tao = PP_inv*B0*B0'*PP_inv;


%%

 L = [ 4 -1 -1 -1  0 -1;
      -1  2 -1  0  0  0;
      -1 -1  3  0 -1  0;
      -1  0  0  2  0 -1;
       0  0 -1  0  2 -1;
      -1  0  0 -1 -1  3];
tmax=200; % maximum simulation time s
T=0.01; % sample time s

kmax=1+tmax/T;


%% 6 agents
zeta1=zeros(Ord*kmax,1);
zeta2=zeros(Ord*kmax,1);
zeta3=zeros(Ord*kmax,1);
zeta4=zeros(Ord*kmax,1);
zeta5=zeros(Ord*kmax,1);
zeta6=zeros(Ord*kmax,1);

% extended states of the agents
x1 = zeros((Ord+1)*kmax,1);
x2 = zeros((Ord+1)*kmax,1);
x3 = zeros((Ord+1)*kmax,1);
x4 = zeros((Ord+1)*kmax,1);
x5 = zeros((Ord+1)*kmax,1);
x6 = zeros((Ord+1)*kmax,1);

% estimation of the states
zeta_hat1=zeros(Ord*kmax,1);
zeta_hat2=zeros(Ord*kmax,1);
zeta_hat3=zeros(Ord*kmax,1);
zeta_hat4=zeros(Ord*kmax,1);
zeta_hat5=zeros(Ord*kmax,1);
zeta_hat6=zeros(Ord*kmax,1);

% estimation of the extended states
x_hat1=zeros((Ord+1)*kmax,1);
x_hat2=zeros((Ord+1)*kmax,1);
x_hat3=zeros((Ord+1)*kmax,1);
x_hat4=zeros((Ord+1)*kmax,1);
x_hat5=zeros((Ord+1)*kmax,1);
x_hat6=zeros((Ord+1)*kmax,1);

% nonlinear part
f1 = zeros(kmax,1);
f2 = zeros(kmax,1);
f3 = zeros(kmax,1);
f4 = zeros(kmax,1);
f5 = zeros(kmax,1);
f6 = zeros(kmax,1);

% disturbance
w1 = zeros(kmax,1);
w2 = zeros(kmax,1);
w3 = zeros(kmax,1);
w4 = zeros(kmax,1);
w5 = zeros(kmax,1);
w6 = zeros(kmax,1);

% coupling weight for agent
num=0;
d1=num*ones(kmax,1);
d2=num*ones(kmax,1);
d3=num*ones(kmax,1);
d4=num*ones(kmax,1);
d5=num*ones(kmax,1);
d6=num*ones(kmax,1);

%% 开始计算，第一次迭代运算
    zeta1(1:Ord)=[5.2873;0.1];
    zeta2(1:Ord)=[0.8778;0.1];
    zeta3(1:Ord)=[1.6757;0.1];
    zeta4(1:Ord)=[1.8980;0.1];
    zeta5(1:Ord)=[3.3125;0.1];
    zeta6(1:Ord)=[0.3381;0.1];
    
%     zeta_hat1(1:Ord)=[5.2873;0.1];
%     zeta_hat2(1:Ord)=[0.8778;0.1];
%     zeta_hat3(1:Ord)=[1.6757;0.1];
%     zeta_hat4(1:Ord)=[1.8980;0.1];
%     zeta_hat5(1:Ord)=[3.3125;0.1];
%     zeta_hat6(1:Ord)=[0.3381;0.1];
    
% zeta_hat1(1:Ord)=zeta1(1:Ord);
% zeta_hat2(1:Ord)=zeta2(1:Ord);
% zeta_hat3(1:Ord)=zeta3(1:Ord);
% zeta_hat4(1:Ord)=zeta4(1:Ord);
% zeta_hat5(1:Ord)=zeta5(1:Ord);
% zeta_hat6(1:Ord)=zeta6(1:Ord);

    
w1(1:1) = 0.3*sin(zeta1(1,1));
w2(1:1) = 0.3*sin(zeta2(1,1));
w3(1:1) = 0.3*sin(zeta3(1,1));
w4(1:1) = 0.3*sin(zeta4(1,1));
w5(1:1) = 0.3*sin(zeta5(1,1));
w6(1:1) = 0.3*sin(zeta6(1,1));

f1(1:1) = sin(zeta1(1,1)) + w1(1:1);
f2(1:1) = sin(zeta2(1,1)) + w2(1:1);
f3(1:1) = sin(zeta3(1,1)) + w3(1:1);
f4(1:1) = sin(zeta4(1,1)) + w4(1:1);
f5(1:1) = sin(zeta5(1,1)) + w5(1:1);
f6(1:1) = sin(zeta6(1,1)) + w6(1:1);

%%
h=waitbar(0,'Please wait...');
for k=1:kmax-1
    waitbar(k/kmax,h);
    
%     w1(k+1) = w1(k) + T*0.3*zeta1(k*Ord)*cos(zeta1(k*(Ord-1)));
%     w2(k+1) = w2(k) + T*0.3*zeta2(k*Ord)*cos(zeta2(k*(Ord-1)));
%     w3(k+1) = w3(k) + T*0.3*zeta3(k*Ord)*cos(zeta3(k*(Ord-1)));
%     w4(k+1) = w4(k) + T*0.3*zeta4(k*Ord)*cos(zeta4(k*(Ord-1)));
%     w5(k+1) = w5(k) + T*0.3*zeta5(k*Ord)*cos(zeta5(k*(Ord-1)));
%     w6(k+1) = w6(k) + T*0.3*zeta6(k*Ord)*cos(zeta6(k*(Ord-1)));
    
    f1(k+1) = f1(k) + T*(1.3*zeta1(k*Ord)*cos(zeta1(Ord*(k-1)+1)));
    f2(k+1) = f2(k) + T*(1.3*zeta2(k*Ord)*cos(zeta2(Ord*(k-1)+1)));
    f3(k+1) = f3(k) + T*(1.3*zeta3(k*Ord)*cos(zeta3(Ord*(k-1)+1)));
    f4(k+1) = f4(k) + T*(1.3*zeta4(k*Ord)*cos(zeta4(Ord*(k-1)+1)));
    f5(k+1) = f5(k) + T*(1.3*zeta5(k*Ord)*cos(zeta5(Ord*(k-1)+1)));
    f6(k+1) = f6(k) + T*(1.3*zeta6(k*Ord)*cos(zeta6(Ord*(k-1)+1)));
    
    inter1 = (L(1,1)*zeta_hat1(Ord*(k-1)+1:Ord*k)+L(1,2)*zeta_hat2(Ord*(k-1)+1:Ord*k)+L(1,3)*zeta_hat3(Ord*(k-1)+1:Ord*k)+L(1,4)*zeta_hat4(Ord*(k-1)+1:Ord*k)+L(1,5)*zeta_hat5(Ord*(k-1)+1:Ord*k)+L(1,6)*zeta_hat6(Ord*(k-1)+1:Ord*k));
    inter2 = (L(2,1)*zeta_hat1(Ord*(k-1)+1:Ord*k)+L(2,2)*zeta_hat2(Ord*(k-1)+1:Ord*k)+L(2,3)*zeta_hat3(Ord*(k-1)+1:Ord*k)+L(2,4)*zeta_hat4(Ord*(k-1)+1:Ord*k)+L(2,5)*zeta_hat5(Ord*(k-1)+1:Ord*k)+L(2,6)*zeta_hat6(Ord*(k-1)+1:Ord*k));
    inter3 = (L(3,1)*zeta_hat1(Ord*(k-1)+1:Ord*k)+L(3,2)*zeta_hat2(Ord*(k-1)+1:Ord*k)+L(3,3)*zeta_hat3(Ord*(k-1)+1:Ord*k)+L(3,4)*zeta_hat4(Ord*(k-1)+1:Ord*k)+L(3,5)*zeta_hat5(Ord*(k-1)+1:Ord*k)+L(3,6)*zeta_hat6(Ord*(k-1)+1:Ord*k));
    inter4 = (L(4,1)*zeta_hat1(Ord*(k-1)+1:Ord*k)+L(4,2)*zeta_hat2(Ord*(k-1)+1:Ord*k)+L(4,3)*zeta_hat3(Ord*(k-1)+1:Ord*k)+L(4,4)*zeta_hat4(Ord*(k-1)+1:Ord*k)+L(4,5)*zeta_hat5(Ord*(k-1)+1:Ord*k)+L(4,6)*zeta_hat6(Ord*(k-1)+1:Ord*k));
    inter5 = (L(5,1)*zeta_hat1(Ord*(k-1)+1:Ord*k)+L(5,2)*zeta_hat2(Ord*(k-1)+1:Ord*k)+L(5,3)*zeta_hat3(Ord*(k-1)+1:Ord*k)+L(5,4)*zeta_hat4(Ord*(k-1)+1:Ord*k)+L(5,5)*zeta_hat5(Ord*(k-1)+1:Ord*k)+L(5,6)*zeta_hat6(Ord*(k-1)+1:Ord*k));
    inter6 = (L(6,1)*zeta_hat1(Ord*(k-1)+1:Ord*k)+L(6,2)*zeta_hat2(Ord*(k-1)+1:Ord*k)+L(6,3)*zeta_hat3(Ord*(k-1)+1:Ord*k)+L(6,4)*zeta_hat4(Ord*(k-1)+1:Ord*k)+L(6,5)*zeta_hat5(Ord*(k-1)+1:Ord*k)+L(6,6)*zeta_hat6(Ord*(k-1)+1:Ord*k));
    
%     u1 = d1(k) * K * (L(1,1)*zeta_hat1(Ord*(k-1)+1:Ord*k)+L(1,2)*zeta_hat2(Ord*(k-1)+1:Ord*k)+L(1,3)*zeta_hat3(Ord*(k-1)+1:Ord*k)+L(1,4)*zeta_hat4(Ord*(k-1)+1:Ord*k)+L(1,5)*zeta_hat5(Ord*(k-1)+1:Ord*k)+L(1,6)*zeta_hat6(Ord*(k-1)+1:Ord*k)) - x_hat1((Ord+1)*k);
%     u2 = d2(k) * K * (L(2,1)*zeta_hat1(Ord*(k-1)+1:Ord*k)+L(2,2)*zeta_hat2(Ord*(k-1)+1:Ord*k)+L(2,3)*zeta_hat3(Ord*(k-1)+1:Ord*k)+L(2,4)*zeta_hat4(Ord*(k-1)+1:Ord*k)+L(2,5)*zeta_hat5(Ord*(k-1)+1:Ord*k)+L(2,6)*zeta_hat6(Ord*(k-1)+1:Ord*k)) - x_hat2((Ord+1)*k);
%     u3 = d3(k) * K * (L(3,1)*zeta_hat1(Ord*(k-1)+1:Ord*k)+L(3,2)*zeta_hat2(Ord*(k-1)+1:Ord*k)+L(3,3)*zeta_hat3(Ord*(k-1)+1:Ord*k)+L(3,4)*zeta_hat4(Ord*(k-1)+1:Ord*k)+L(3,5)*zeta_hat5(Ord*(k-1)+1:Ord*k)+L(3,6)*zeta_hat6(Ord*(k-1)+1:Ord*k)) - x_hat3((Ord+1)*k);
%     u4 = d4(k) * K * (L(4,1)*zeta_hat1(Ord*(k-1)+1:Ord*k)+L(4,2)*zeta_hat2(Ord*(k-1)+1:Ord*k)+L(4,3)*zeta_hat3(Ord*(k-1)+1:Ord*k)+L(4,4)*zeta_hat4(Ord*(k-1)+1:Ord*k)+L(4,5)*zeta_hat5(Ord*(k-1)+1:Ord*k)+L(4,6)*zeta_hat6(Ord*(k-1)+1:Ord*k)) - x_hat4((Ord+1)*k);
%     u5 = d5(k) * K * (L(5,1)*zeta_hat1(Ord*(k-1)+1:Ord*k)+L(5,2)*zeta_hat2(Ord*(k-1)+1:Ord*k)+L(5,3)*zeta_hat3(Ord*(k-1)+1:Ord*k)+L(5,4)*zeta_hat4(Ord*(k-1)+1:Ord*k)+L(5,5)*zeta_hat5(Ord*(k-1)+1:Ord*k)+L(5,6)*zeta_hat6(Ord*(k-1)+1:Ord*k)) - x_hat5((Ord+1)*k);
%     u6 = d6(k) * K * (L(6,1)*zeta_hat1(Ord*(k-1)+1:Ord*k)+L(6,2)*zeta_hat2(Ord*(k-1)+1:Ord*k)+L(6,3)*zeta_hat3(Ord*(k-1)+1:Ord*k)+L(6,4)*zeta_hat4(Ord*(k-1)+1:Ord*k)+L(6,5)*zeta_hat5(Ord*(k-1)+1:Ord*k)+L(6,6)*zeta_hat6(Ord*(k-1)+1:Ord*k)) - x_hat6((Ord+1)*k);

    u1 = d1(k) * K * inter1 - x_hat1((Ord+1)*k);
    u2 = d2(k) * K * inter2 - x_hat2((Ord+1)*k);
    u3 = d3(k) * K * inter3 - x_hat3((Ord+1)*k);
    u4 = d4(k) * K * inter4 - x_hat4((Ord+1)*k);
    u5 = d5(k) * K * inter5 - x_hat5((Ord+1)*k);
    u6 = d6(k) * K * inter6 - x_hat6((Ord+1)*k);
   
    zeta1(Ord*k+1:Ord*(k+1)) = zeta1(Ord*(k-1)+1:Ord*k) + T*(A0*zeta1(Ord*(k-1)+1:Ord*k)+B0*(u1 + f1(k)));
    zeta2(Ord*k+1:Ord*(k+1)) = zeta2(Ord*(k-1)+1:Ord*k) + T*(A0*zeta2(Ord*(k-1)+1:Ord*k)+B0*(u2 + f2(k)));
    zeta3(Ord*k+1:Ord*(k+1)) = zeta3(Ord*(k-1)+1:Ord*k) + T*(A0*zeta3(Ord*(k-1)+1:Ord*k)+B0*(u3 + f3(k)));
    zeta4(Ord*k+1:Ord*(k+1)) = zeta4(Ord*(k-1)+1:Ord*k) + T*(A0*zeta4(Ord*(k-1)+1:Ord*k)+B0*(u4 + f4(k)));
    zeta5(Ord*k+1:Ord*(k+1)) = zeta5(Ord*(k-1)+1:Ord*k) + T*(A0*zeta5(Ord*(k-1)+1:Ord*k)+B0*(u5 + f5(k)));
    zeta6(Ord*k+1:Ord*(k+1)) = zeta6(Ord*(k-1)+1:Ord*k) + T*(A0*zeta6(Ord*(k-1)+1:Ord*k)+B0*(u6 + f6(k)));
    

      zeta_hat1(Ord*k+1:Ord*(k+1)) = zeta_hat1(Ord*(k-1)+1:Ord*k) + T*(A0*zeta_hat1(Ord*(k-1)+1:Ord*k)+B0*(u1 + x_hat1((Ord+1)*k))+G0_star*C0*(zeta1(Ord*(k-1)+1:Ord*k)-zeta_hat1(Ord*(k-1)+1:Ord*k)));
      zeta_hat2(Ord*k+1:Ord*(k+1)) = zeta_hat2(Ord*(k-1)+1:Ord*k) + T*(A0*zeta_hat2(Ord*(k-1)+1:Ord*k)+B0*(u2 + x_hat2((Ord+1)*k))+G0_star*C0*(zeta2(Ord*(k-1)+1:Ord*k)-zeta_hat2(Ord*(k-1)+1:Ord*k)));
      zeta_hat3(Ord*k+1:Ord*(k+1)) = zeta_hat3(Ord*(k-1)+1:Ord*k) + T*(A0*zeta_hat3(Ord*(k-1)+1:Ord*k)+B0*(u3 + x_hat3((Ord+1)*k))+G0_star*C0*(zeta3(Ord*(k-1)+1:Ord*k)-zeta_hat3(Ord*(k-1)+1:Ord*k)));
      zeta_hat4(Ord*k+1:Ord*(k+1)) = zeta_hat4(Ord*(k-1)+1:Ord*k) + T*(A0*zeta_hat4(Ord*(k-1)+1:Ord*k)+B0*(u4 + x_hat4((Ord+1)*k))+G0_star*C0*(zeta4(Ord*(k-1)+1:Ord*k)-zeta_hat4(Ord*(k-1)+1:Ord*k)));
      zeta_hat5(Ord*k+1:Ord*(k+1)) = zeta_hat5(Ord*(k-1)+1:Ord*k) + T*(A0*zeta_hat5(Ord*(k-1)+1:Ord*k)+B0*(u5 + x_hat5((Ord+1)*k))+G0_star*C0*(zeta5(Ord*(k-1)+1:Ord*k)-zeta_hat5(Ord*(k-1)+1:Ord*k)));
      zeta_hat6(Ord*k+1:Ord*(k+1)) = zeta_hat6(Ord*(k-1)+1:Ord*k) + T*(A0*zeta_hat6(Ord*(k-1)+1:Ord*k)+B0*(u6 + x_hat6((Ord+1)*k))+G0_star*C0*(zeta6(Ord*(k-1)+1:Ord*k)-zeta_hat6(Ord*(k-1)+1:Ord*k)));
      
      x_hat1((Ord+1)*(k+1)) = x_hat1((Ord+1)*k) + T*(w0^3*alpha3*(zeta1(Ord*(k-1)+1)-zeta_hat1(Ord*(k-1)+1))+1.3*zeta_hat1(k*Ord)*cos(zeta_hat1(k*(Ord-1))));
      x_hat2((Ord+1)*(k+1)) = x_hat2((Ord+1)*k) + T*(w0^3*alpha3*(zeta2(Ord*(k-1)+1)-zeta_hat2(Ord*(k-1)+1))+1.3*zeta_hat2(k*Ord)*cos(zeta_hat2(k*(Ord-1))));
      x_hat3((Ord+1)*(k+1)) = x_hat3((Ord+1)*k) + T*(w0^3*alpha3*(zeta3(Ord*(k-1)+1)-zeta_hat3(Ord*(k-1)+1))+1.3*zeta_hat3(k*Ord)*cos(zeta_hat3(k*(Ord-1))));
      x_hat4((Ord+1)*(k+1)) = x_hat4((Ord+1)*k) + T*(w0^3*alpha3*(zeta4(Ord*(k-1)+1)-zeta_hat4(Ord*(k-1)+1))+1.3*zeta_hat4(k*Ord)*cos(zeta_hat4(k*(Ord-1))));
      x_hat5((Ord+1)*(k+1)) = x_hat5((Ord+1)*k) + T*(w0^3*alpha3*(zeta5(Ord*(k-1)+1)-zeta_hat5(Ord*(k-1)+1))+1.3*zeta_hat5(k*Ord)*cos(zeta_hat5(k*(Ord-1))));
      x_hat6((Ord+1)*(k+1)) = x_hat6((Ord+1)*k) + T*(w0^3*alpha3*(zeta6(Ord*(k-1)+1)-zeta_hat6(Ord*(k-1)+1))+1.3*zeta_hat6(k*Ord)*cos(zeta_hat6(k*(Ord-1))));
      
      dotd1 = kappa*inter1'*tao*inter1;
      dotd2 = kappa*inter2'*tao*inter2;
      dotd3 = kappa*inter3'*tao*inter3;
      dotd4 = kappa*inter4'*tao*inter4;
      dotd5 = kappa*inter5'*tao*inter5;
      dotd6 = kappa*inter6'*tao*inter6;
      
      d1(k+1) = d1(k) + T*dotd1;
      d2(k+1) = d2(k) + T*dotd2;
      d3(k+1) = d3(k) + T*dotd3;
      d4(k+1) = d4(k) + T*dotd4;
      d5(k+1) = d5(k) + T*dotd5;
      d6(k+1) = d6(k) + T*dotd6;
      
end

close(h);


figure(1);
plot(0:T:T*(kmax-1),zeta2(1:Ord:length(zeta2))-zeta1(1:Ord:length(zeta1)));
hold on;
plot(0:T:T*(kmax-1),zeta3(1:Ord:length(zeta3))-zeta1(1:Ord:length(zeta1)));
hold on;
plot(0:T:T*(kmax-1),zeta4(1:Ord:length(zeta4))-zeta1(1:Ord:length(zeta1)));
hold on;
plot(0:T:T*(kmax-1),zeta5(1:Ord:length(zeta5))-zeta1(1:Ord:length(zeta1)));
hold on;
plot(0:T:T*(kmax-1),zeta6(1:Ord:length(zeta6))-zeta1(1:Ord:length(zeta1)));
hold on;
plot(0:T:T*(kmax-1),zeta2(2:Ord:length(zeta2))-zeta1(2:Ord:length(zeta1)));
hold on;
plot(0:T:T*(kmax-1),zeta3(2:Ord:length(zeta3))-zeta1(2:Ord:length(zeta1)));
hold on;
plot(0:T:T*(kmax-1),zeta4(2:Ord:length(zeta4))-zeta1(2:Ord:length(zeta1)));
hold on;
plot(0:T:T*(kmax-1),zeta5(2:Ord:length(zeta5))-zeta1(2:Ord:length(zeta1)));
hold on;
plot(0:T:T*(kmax-1),zeta6(2:Ord:length(zeta6))-zeta1(2:Ord:length(zeta1)));
hold on;
grid on;
xlabel('time $(s)$','Interpreter','latex');
ylabel('$x_{i}-x_{1}$','Interpreter','latex');
print -depsc -r600 errors

figure(2)
plot(0:T:T*(kmax-1),zeta_hat1(1:Ord:length(zeta1)));
hold on;
plot(0:T:T*(kmax-1),zeta1(1:Ord:length(zeta1)));
hold on;
grid on;

figure(3)
plot(0:T:T*(kmax-1),zeta_hat1(1:Ord:length(zeta_hat1))-zeta1(1:Ord:length(zeta1)));
grid on;
plot(0:T:T*(kmax-1),zeta_hat1(2:Ord:length(zeta_hat1))-zeta1(2:Ord:length(zeta1)));
grid on;
hold on;

figure(4)
plot(0:T:T*(kmax-1),zeta_hat2(1:Ord:length(zeta_hat2))-zeta2(1:Ord:length(zeta2)));
grid on;
plot(0:T:T*(kmax-1),zeta_hat2(2:Ord:length(zeta_hat2))-zeta2(2:Ord:length(zeta2)));
grid on;
hold on;

figure(5)
plot(0:T:T*(kmax-1),zeta_hat3(1:Ord:length(zeta_hat3))-zeta3(1:Ord:length(zeta3)));
grid on;
hold on;

figure(6)
plot(0:T:T*(kmax-1),zeta_hat4(1:Ord:length(zeta_hat4))-zeta4(1:Ord:length(zeta4)));
grid on;
hold on;

figure(7)
plot(0:T:T*(kmax-1),zeta_hat5(1:Ord:length(zeta_hat5))-zeta5(1:Ord:length(zeta5)));
grid on;
hold on;

figure(8)
plot(0:T:T*(kmax-1),zeta_hat6(1:Ord:length(zeta_hat6))-zeta6(1:Ord:length(zeta6)));
grid on;
hold on;

figure(9)
subplot(3,2,1)
plot(0:T:T*(kmax-1),zeta_hat1(1:Ord:length(zeta_hat1))-zeta1(1:Ord:length(zeta1)));
grid on;
plot(0:T:T*(kmax-1),zeta_hat1(2:Ord:length(zeta_hat1))-zeta1(2:Ord:length(zeta1)));
grid on;
hold on;

subplot(3,2,2)
plot(0:T:T*(kmax-1),zeta_hat2(1:Ord:length(zeta_hat2))-zeta2(1:Ord:length(zeta2)));
grid on;
plot(0:T:T*(kmax-1),zeta_hat2(2:Ord:length(zeta_hat2))-zeta2(2:Ord:length(zeta2)));
grid on;
hold on;

subplot(3,2,3)
plot(0:T:T*(kmax-1),zeta_hat3(1:Ord:length(zeta_hat3))-zeta3(1:Ord:length(zeta3)));
grid on;
hold on;
plot(0:T:T*(kmax-1),zeta_hat3(2:Ord:length(zeta_hat3))-zeta3(2:Ord:length(zeta3)));
hold on;

subplot(3,2,4)
plot(0:T:T*(kmax-1),zeta_hat4(1:Ord:length(zeta_hat4))-zeta4(1:Ord:length(zeta4)));
hold on;
plot(0:T:T*(kmax-1),zeta_hat4(2:Ord:length(zeta_hat4))-zeta4(2:Ord:length(zeta4)));
hold on;
grid on;

subplot(3,2,5)
plot(0:T:T*(kmax-1),zeta_hat5(1:Ord:length(zeta_hat5))-zeta5(1:Ord:length(zeta5)));
hold on;
plot(0:T:T*(kmax-1),zeta_hat5(2:Ord:length(zeta_hat5))-zeta5(2:Ord:length(zeta5)));
hold on;
grid on;

subplot(3,2,6)
plot(0:T:T*(kmax-1),zeta_hat6(1:Ord:length(zeta_hat6))-zeta6(1:Ord:length(zeta6)));
hold on;
plot(0:T:T*(kmax-1),zeta_hat6(2:Ord:length(zeta_hat6))-zeta6(2:Ord:length(zeta6)));
hold on;
grid on;
hold on;



        