function pos=estimate_pos(Y,init_pos)

[M,N]=size(Y);

Ys=Y(1:M/2,:);
s=[mean(Ys(1:3:end,:)); mean(Ys(2:3:end,:)); mean(Ys(3:3:end,:))];
Ys=Ys-repmat(s,M/6,1);

Yw=Y(M/2+1:end,:);
w=[mean(Yw(1:3:end,:)); mean(Yw(2:3:end,:)); mean(Yw(3:3:end,:))];
wacc = (w(:,3:end) - w(:,1:end-2))/(2*0.001); 
wacc = [zeros(3,1) wacc zeros(3,1)];



x=reshape(init_pos,M/2,1);
P=0.001^2*eye(M/2);
R=0.01^2*eye(M/2);


for n=2:N-1
    if all(abs(Yw(:,n-1:n+1))>1  & abs(Yw(:,n-1:n+1))<1800*pi/180)
    H=kron(eye(M/6),(skew_sym(w(:,n))^2+skew_sym(wacc(:,n))));
    Re=(H*P*H'+R);
    K=(P*H')/Re;
    x=x+K*(Ys(:,n)-H*x);
    P=(eye(M/2)-K*H)*P;
    end
end
pos=reshape(x,size(init_pos));
end


