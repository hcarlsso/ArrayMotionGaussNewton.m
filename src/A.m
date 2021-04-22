function z=A(u,v)
z=skew_sym(u)'*skew_sym(v)+skew_sym(my_cross(v,u));
end

