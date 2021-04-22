function J=CostFunc2(x,t,y)

N=length(t);

w=x(1).*t+x(2).*sin(x(1).*t+x(3)*ones(N,1));

J=y-x(4)*sin(w.*t+x(5))+x(6)*ones(N,1);


end

