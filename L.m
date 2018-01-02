function [sys,x0,str,ts]=L(t,x,u,flag)

switch flag
  case 0
    [sys,x0,str,ts]=mdlInitializeSizes;
  case 1
    sys=mdlDerivatives(t,x,u);
  case 3
    sys=mdlOutputs(t,x,u);
  case{1,2,4,9}
    sys=[];
  otherwise
    error(['Unhandled flag=',num2str(flag)]);

end

function [sys,x0,str,ts]=mdlInitializeSizes
sizes = simsizes;

sizes.NumContStates  = 3;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 2;
sizes.NumInputs      = 3;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 0;   % at least one sample time is needed

sys = simsizes(sizes);

x0  = [0;0;0];
str = [];
ts  = [];

function sys=mdlDerivatives(t,x,u)
sys(1) = 0;
sys(2) = 0;
sys(3) = 0;


function sys=mdlOutputs(t,x,u)
sys(1) = 40*x(1)-20*x(3);
sys(2) = -20*x(2)+60*x(3);