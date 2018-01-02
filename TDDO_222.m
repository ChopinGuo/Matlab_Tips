function [sys,x0,str,ts]=TDDO_222(t,x,u,flag)

switch flag
  case 0
    [sys,x0,str,ts]=mdlInitializeSizes;
  case 1
    sys=mdlDerivatives(t,x,u);
  case 3
    sys=mdlOutputs(t,x,u);
  case{2,4,9}
    sys=[];
  otherwise
    error(['Unhandled flag=',num2str(flag)]);

end

function [sys,x0,str,ts]=mdlInitializeSizes
sizes = simsizes;

sizes.NumContStates  = 2;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 2;
sizes.NumInputs      = 5;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 0;   % at least one sample time is needed

sys = simsizes(sizes);

x0  = [0;0];
str = [];
ts  = [];

function sys=mdlDerivatives(t,x,u)

a1 = u(1);
a2 = u(2);
b1 = u(3);
b2 = u(4);
b3 = u(5);
% z1 = x(1);
% z2 = x(2);

sys(1) = 20*a2 - 1568*b1 - 520*b2 + 2056*b3 - 40*x(1) + 24*x(2);
sys(2) = 20*a1 - 60*a2 + 960*b1 + 1780*b2 - 5600*b3 + 24*x(1) - 92*x(2);



function sys=mdlOutputs(t,x,u)
sys(1) = x(1);
sys(2) = x(2);
