function[sys,x0]=s_function(t,x,u,flag)
if flag==0
    sys=[2,0,1,1,0,0];
    x0=[0,0];
elseif flag==1
    sys(1)=x(2);
    sys(2)=-(25+10*rand(1))*x(2)+(133+30*rand(1))*u;
elseif flag==3
    sys(1)=x(1);
else
    sys=[];
end