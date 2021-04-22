function J=Grad(w,Settings)
J=repmat(eye(3),Settings.Ns+Settings.Nw,1); 

for k=1:Settings.Ns
    J(3*(k-1)+1:3*k,1:3)=A(w,Settings.R(:,k));
end
end

function z=A(u,v)
z=skew_sym(u)'*skew_sym(v)+skew_sym(my_cross(v,u));
end

function z=my_cross(u,v)
z=[u(2)*v(3)-u(3)*v(2); u(3)*v(1)-u(1)*v(3);u(1)*v(2)-u(2)*v(1)];
end