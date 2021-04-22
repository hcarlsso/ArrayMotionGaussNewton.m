function CRB=CRB_general_case(Settings,w)

% Allocate memory
Settings.N=size(w,2);
CRB=zeros(9,Settings.N);


% Generate accelerometer location matrix
Settings.R=Settings.alpha*Settings.acc_geometry;

% Generate the H matrix in the signal model.
H=zeros(3*(Settings.Ns+Settings.Nw),6);
for k=1:Settings.Ns
   H(3*(k-1)+1:3*k,1:3)=-skew_sym(Settings.R(:,k));
   H(3*(k-1)+1:3*k,4:6)=eye(3);  
end

% Generate the measurement error covariance
% Generate the measurement error
Qs=Settings.sigma_s^2*eye(3*Settings.Ns);
Qw=(Settings.sigma_w*pi/180)^2*eye(3*Settings.Nw);
Q=blkdiag(Qs,Qw);
Qinv=inv(Q);

for n=1:Settings.N
    
    % Get thhe current angular velocity vector
    w_rad=pi/180*w(:,n)*Settings.w_norm;
 
    
    % Check if any gyroscopes are saturated and generate the associated
    % indicator matrix to remove these entries from the signal model
    ind=logical([ones(3*Settings.Ns,1); ...
        repmat((w_rad<=pi/180*Settings.gamma_w),Settings.Nw,1)]);
    
    % Calculate the gradient of h(w) and remove the entries of the 
    % saturated gyroscopes
    J=Grad(w_rad,Settings);
    
    % Calculate the Fisher information matrix
    Phi=[J H]; 
    I=Phi(ind,:)'*Qinv(ind,ind)*Phi(ind,:);
    
    % Invert the Fisher information matrix
    CRB(:,n)=diag(inv(I));
end


end

