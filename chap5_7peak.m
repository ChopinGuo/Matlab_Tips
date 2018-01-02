function [sys,x0,str,ts] = chap5_7peak(t,x,u,flag)
switch flag
  case 0
    [sys,x0,str,ts]=mdlInitializeSizes;
  case 3
    sys=mdlOutputs(t,x,u);
  case{1,2,4,9}
    sys=[];
  otherwise
    DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));

end

function [sys,x0,str,ts]=mdlInitializeSizes

sizes = simsizes;

sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 2;
sizes.NumInputs      = 0;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);

x0  = [];
str = [];
ts  = [0 0];

function sys=mdlOutputs(t,x,u)
Lamda = 50;
R = 100*(1-exp(-Lamda*t)) / (1+exp(-Lamda*t));
Epsilon = 1/R;
sys(1) = R;
sys(2) = Epsilon;
