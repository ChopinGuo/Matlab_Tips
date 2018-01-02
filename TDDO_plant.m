function [sys,x0,str,ts]=TDDO_plant(t,x,u,flag)

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

function [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes
sizes = simsizes;

sizes.NumContStates  = 3;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 5;
sizes.NumInputs      = 2;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 0;   % at least one sample time is needed

sys = simsizes(sizes);

x0  = [0;0;0];
str = [];
ts  = [];

function sys=mdlDerivatives(t,x,u)

% d1=0*sin(0.001*t);
% d2=0*sin(0.001*t);

d1 = 2;
d2 = 3;

sys(1) = -0.8*x(1) + x(2) + 1.6*x(3) + 0.8*d1;
sys(2) = -3*x(2) + 2*x(3) + u(1) - d2;
sys(3) = -6*x(3) + u(2) -0.4*d1 + 1.2*d2;


function sys=mdlOutputs(t,x,u)
% d1=150*sign(sin(0.1*t));
% d2=120*sign(sin(0.1*t));
d1 = 2;
d2 = 3;
sys(1) = x(1);
sys(2) = x(2);
sys(3) = x(3);
sys(4) = d1;
sys(5) = d2;