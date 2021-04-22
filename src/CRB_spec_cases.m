function [CRB1,CRB2,CRB_lim]=CRB_spec_cases(Settings,w)

w_rad=w;

% Allocate memory
CRB1=zeros(9,Settings.N);
CRB2=zeros(9,Settings.N);
CRB_lim=zeros(3,2);

% Calculate the normalized Gamma matrices
Gamma_bar_11=zeros(3);
Gamma_bar_12=zeros(3);
Gamma_bar_22=zeros(3);

for n=1:Settings.Ns
Gamma_bar_11=Gamma_bar_11+A(Settings.w_norm,Settings.acc_geometry(:,n))'...
    *A(Settings.w_norm,Settings.acc_geometry(:,n));

Gamma_bar_12=Gamma_bar_12+A(Settings.w_norm,Settings.acc_geometry(:,n))'...
    *skew_sym(Settings.acc_geometry(:,n));

Gamma_bar_22=Gamma_bar_22+skew_sym(Settings.acc_geometry(:,n))'...
    *skew_sym(Settings.acc_geometry(:,n));
end

Gamma_bar_11=Gamma_bar_11./Settings.Ns;
Gamma_bar_12=-Gamma_bar_12./Settings.Ns;
Gamma_bar_22=Gamma_bar_22./Settings.Ns;

% Schur complements
Schur11=Gamma_bar_22-Gamma_bar_12'*(Gamma_bar_11\Gamma_bar_12);
Schur22=Gamma_bar_11-Gamma_bar_12*(Gamma_bar_22\Gamma_bar_12');



%% CRB for special case 1 
for n=1:Settings.N
    CRB1(1:3,n)=(Settings.sigma_w)^2/Settings.Nw*diag(inv(eye(3)+...
        Settings.alpha^2*w_rad(n)^2*((Settings.sigma_w)^2/Settings.Nw)...
        /(Settings.sigma_s^2/Settings.Ns)*Schur22));
    
    CRB1(4:6,n)=Settings.sigma_s^2/(Settings.alpha^2*Settings.Ns)*diag(inv(...
        Gamma_bar_22)-Settings.alpha^2*w_rad(n)^2*(Settings.sigma_w^2/Settings.Nw)...
        /(Settings.sigma_s^2/Settings.Ns)*Gamma_bar_12'*inv(eye(3)+Settings.alpha^2*w_rad(n)^2*(Settings.sigma_w^2/Settings.Nw)...
        /(Settings.sigma_s^2/Settings.Ns)*Gamma_bar_11)*Gamma_bar_12);
end  

CRB1(7:9,:)=Settings.sigma_s^2/Settings.Ns*ones(3,Settings.N);

% Limit values for the CRB for the angular acceleration
CRB_lim(:,1)=Settings.sigma_s^2/(Settings.alpha^2*Settings.Ns)*diag(inv(...
    Gamma_bar_22));

CRB_lim(:,2)=Settings.sigma_s^2/(Settings.alpha^2*Settings.Ns)*diag(inv(...
    Schur11));


%% CRB for special case 2
CRB2(1:3,:)=Settings.sigma_s^2/(Settings.Ns*Settings.alpha^2)*...
    diag(inv(Schur22))*(1./(w_rad.^2));

CRB2(4:6,:)=Settings.sigma_s^2/(Settings.Ns*Settings.alpha^2)*...
    diag(inv(Schur11))*ones(1,Settings.N);

CRB2(7:9,:)=Settings.sigma_s^2/Settings.Ns*ones(3,Settings.N);


end


function z=A(u,v)
z=skew_sym(u)'*skew_sym(v)+skew_sym(cross(v,u));
end
