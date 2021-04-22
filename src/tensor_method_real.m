function x=tensor_method_real(y,Settings)


% Allocate memory
x=zeros(9,1);

% Generate accelerometer location matrix
R=[Settings.alpha*Settings.acc_geometry; ones(1,Settings.Ns)];

% Generate projection matrix
P=R'/(R*R');

% Least square estimate
T=reshape(y(1:3*Settings.Ns),3,Settings.Ns)*P;

% Angular velocity
Ws=0.5*(T(1:3,1:3)+T(1:3,1:3)');

x(1:3)=sqrt(abs(diag(Ws-0.5*trace(Ws)*eye(3))));

% Angular acceleration
x(4)=0.5*(T(3,2)-T(2,3));
x(5)=0.5*(T(1,3)-T(3,1));
x(6)=0.5*(T(2,1)-T(1,2));

% Specific force
x(7:9)=T(:,4);


end
