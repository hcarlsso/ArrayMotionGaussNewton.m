function [Y,w]=generate_data(Settings)
% Keep the angular accerleration and specific force constant 
% While increasing the angular velocity
    
% Angular speed vector (degrees,radians,and truncated radians)
%w=linspace(Settings.wmin,Settings.wmax,Settings.N);
%w=logspace(log10(Settings.wmin),log10(Settings.wmax),Settings.N-2);
a=1/max(Settings.w_norm);
w=[logspace(... % From wmin to max of normlized angular velocity
    log10(Settings.wmin),...
    log10(a*Settings.gamma_w-3),...
    Settings.N/2)...
   logspace(...
       log10(a*Settings.gamma_w+3),...
       log10(Settings.wmax),...
       Settings.N/2)];
w_rad=w.*pi/180;


% Generate the measurement error
Qs=Settings.sigma_s^2*eye(3*Settings.Ns);
Qw=(Settings.sigma_w*pi/180)^2*eye(3*Settings.Nw);
% Construct block diagonal matrix from input arguments
Q=blkdiag(Qs,Qw);

% Cholesky factorization
% R = chol(A) produces an upper triangular matrix R from 
% the diagonal and upper triangle of matrix A, satisfying 
% the equation R'*R=A. The chol function assumes that A is 
% (complex Hermitian) symmetric. If it is not, chol uses the 
% (complex conjugate) transpose of the upper triangle as the 
% lower triangle. Matrix A must be positive definite.

% To get the standard deviation?
N=chol(Q)*randn(3*(Settings.Ns+Settings.Nw),Settings.N);

% Generate accelerometer location matrix
R= Settings.alpha*... % Scalar
   Settings.acc_geometry+... % N_sensors x 3 , size of R
   sqrt(Settings.s2r)*... % Standard deviation of the error in
   ...                       % sensor location, why sqrt of std dev?
   randn(size(Settings.acc_geometry)); % Add some noise to positions

% Generate the h(w) part of the signal model for the normalized angular 
% velocity vector. 
hs_norm=reshape(...
    (...
        skew_sym(...
            Settings.w_norm ... % w_norm : 3x1
            )^2 ... % 3x3
        )*R, ... % 3 x N_sensors , size of R
    3*Settings.Ns,1); % why to make a col vector of it?
                    
hw_norm=repmat(Settings.w_norm,Settings.Nw,1); % Repeat the ang Vec
                                               % Nw times

% Generate the phi vector
phi=[Settings.Angular_Acceleration;Settings.Specific_Force];

% Generate the H matrix in the signal model.
H=zeros(3*(Settings.Ns+Settings.Nw),6);
% Stuff that is not filled is for the gyroscopes
for k=1:Settings.Ns
   H(3*(k-1)+1:3*k,1:3)=-skew(R(:,k)); % Put in Omega_{r_i}
   H(3*(k-1)+1:3*k,4:6)=eye(3);  
end

% Generate the error free signal. 
X=[hs_norm*(w_rad.^2); hw_norm*w_rad]+repmat(H*phi,1,Settings.N);

% Add the measurement noise
Y=X+N; 

% Check if any of the gyroscopes are saturated
TMP=Y(3*Settings.Ns+1:end,:);
TMP(TMP>(Settings.gamma_w*pi/180))=Settings.gamma_w*pi/180;
Y(3*Settings.Ns+1:end,:)=TMP;
end

