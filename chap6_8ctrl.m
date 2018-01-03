function [sys,x0,str,ts] = sfuntmpl(t,x,u,flag)

switch flag
  case 0
    [sys,x0,str,ts]=mdlInitializeSizes;
  case 3
    sys=mdlOutputs(t,x,u);
  case {1,2,4,9}
    sys=[];
  otherwise
    DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));

end

function [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes

sizes = simsizes;

sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 1;
sizes.NumInputs      = 5;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 0;  

sys = simsizes(sizes);
x0  = [];
str = [];
ts  = [];

function sys=mdlOutputs(t,x,u)
yd = u(1);
dyd = cos(t);
y = u(2);
dy = u(3);
fp = u(5);

M = 2;

if M==1
    ut = 10*(yd-y)+10*(dyd-dy);
elseif M==2
    ut = 10*(yd-y)+10*(dyd-dy)-1/133*fp;
end

sys(1) = ut;
