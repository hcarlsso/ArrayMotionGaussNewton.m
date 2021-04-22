function J=CostFunc(x,t,y)
w0=x(1);
b1=x(2);
b2=x(3);

w=w0+b1*sin(w0.*t)+b2*cos(w0.*t);
H=[sin(w.*t) cos(w.*t) ones(size(t))];

J=-(y'*H)*((H'*H)\(H'*y));
end

