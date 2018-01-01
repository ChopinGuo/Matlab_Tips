% S-function for continuous state equation
function [sys,x0,str,ts] = s_funtion(t,x,u,flag)

switch flag,
% Initialization 
  case 0,
    [sys,x0,str,ts]=mdlInitializeSizes;


% Derivatives
  case 1,
    sys=mdlDerivatives(t,x,u);

% Outputs 
  case 3,
    sys=mdlOutputs(t,x,u);

% Unhandled flags
  case {2, 4, 9},
    sys=[];


% Unexpected flags 
  otherwise
    DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));

end

% mdlInitializeSizes
function [sys,x0,str,ts]=mdlInitializeSizes

sizes = simsizes;

sizes.NumContStates  = 2;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 1;
sizes.NumInputs      = 1;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 0;   % at least one sample time is needed

sys = simsizes(sizes);

x0  = [0,0];
str = [];
ts  = [];

function sys=mdlDerivatives(t,x,u)
sys(1) = x(2);
sys(2) = -(25+10*rands(1))*x(2)+(133+30*rand(1))*u;

function sys=mdlOutputs(t,x,u)

sys(1) = x(1);
