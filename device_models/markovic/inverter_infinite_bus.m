function inverter_dxdt = inverter_infinite_bus(t,x)
    parameters;
    eg_d = x(1);
    eg_q = x(2);
    ig_d = x(3);
    ig_q = x(4);
    is_d = x(5);
    is_q = x(6);
    xi_d = x(7);
    xi_q = x(8);
    gamma_d = x(9);
    gamma_q = x(10);
    epsilon = x(11);
    theta_apc = x(12);
    theta_pll= x(13);
    pc_tilde = x(14);
    qc_tilde = x(15);
    pc=x(16);
    qc=x(17);
    eg_hat_d=x(18);
    eg_hat_q=x(19);
    w_pll=x(20);
    wc=x(21);
    vc=x(22);
    vc_bar_d=x(23);
    vc_bar_q=x(24);
    is_bar_d=x(25);
    is_bar_q=x(26);
    vm_d=x(27);
    vm_q=x(28);
        
% Differential Equations
% Eq 5 d(eg_d)/dt
f=[(wb/cf)*(is_d-ig_d)+wc*wb*eg_q];
    
% Eq 5 d(eg_q)/dt 
f = [f; (wb/cf)*(is_q-ig_q)-wc*wb*eg_d];
    
% Eq 4 d(ig_d)/dt 
f = [f; (wb/lf)*(eg_d-vn_d)-(rf*wb/lf)*ig_d+wb*wc*ig_q];
      
% Eq 4 d(ig_q)/dt  
f = [f; (wb/lf)*(eg_q-vn_q)-wb*wc*ig_d-(rf*wb/lf)*ig_q];
      
% Eq 3 d(is_d)/dt  
f = [f; (wb/lf)*(vm_d-eg_d)-(rf*wb/lf)*is_d+wb*wc*is_q];
      
% Eq 3 d(is_q)/dt  
f = [f; (wb/lf)*(vm_q-eg_q)-wb*wc*is_d-(rf*wb/lf)*is_q];

% Eq 21 d(xi_d)/dt
f = [f; vc_bar_d - eg_d];
    
% Eq 21 d(xi_q)/dt
f = [f; vc_bar_q - eg_q];
      
% Eq. 23 d(gamma_d)/dt=
f = [f; is_bar_d-is_d];
    
% Eq. 23 d(gamma_q)/dt=
f = [f; is_bar_q-is_q];

% Eq. 10 d(epsilon)/dt = 
f = [f; eg_hat_q];
      
% Eq. 8 d(theta_apc)/dt = 
f = [f; wc*wb];
      
% Eq. 11 d(theta_pll)/dt=
f = [f; w_pll*wb];
      
% Eq 13. d(pc_tilde)/dt = 
f = [f; 2*pi*fz*(pc - pc_tilde)];
      
% Eq. 17 d(qc_tilde)/dt=
f = [f; 2*pi*fz*(qc-qc_tilde)];

% Algebraic Equations
% Eq. 6 Pc calc
g = [-(pc) + eg_d*ig_d+eg_q*ig_q];

% Eq. 6 Qc calc
g = [g; -(qc)+eg_q*ig_d-eg_d*ig_q];

% Eq. 7 eg_hat_d
g = [g; -(eg_hat_d)-eg_d*cos(theta_pll-theta_apc)+eg_d*sin(theta_pll-theta_apc)];

% Eq. 7 eg_hat_q
g = [g; -(eg_hat_q)-eg_d*sin(theta_pll-theta_apc)-eg_d*cos(theta_pll-theta_apc)];
      
% Eq. 9 wpll
g = [g; -(w_pll)+w0+kp_pll*eg_hat_q+ki_pll*epsilon];

% Eq. 12 wc
g = [g; -(wc)+wc_ref+rp*(pc_ref - pc_tilde)];
    
% Eq. 16 vc
g = [g; -(vc)+vc_ref+rq*(qc_ref - qc_tilde)];

% Eq. 18 vc_bar
g =  [g; -(vc_bar_d)+vc-rv*ig_d+wc*lv*ig_q];
    
% Eq 19 vc_bar
g =  [g;  -(vc_bar_q)-rv*ig_q+wc*lv*ig_d];


% Eq 20 is_bar
g = [g; -(is_bar_d)+kp_v*(vc_bar_d-eg_d)+ki_v*xi_d - wc*cf*eg_q + kf_i*ig_d];
    
g = [g; -(is_bar_q)+kp_v*(vc_bar_d-eg_q)+ki_v*xi_q + wc*cf*eg_d + kf_i*ig_q];
    
% Eq 22 vm_d
g = [g; -(vm_d)+kp_i*(is_bar_d-is_d)+ki_i*gamma_d - wc*lf*is_q + kf_v*eg_d];
    
g = [g; -(vm_q)+kp_i*(is_bar_q-is_q)+ki_i*gamma_q + wc*lf*is_d + kf_v*eg_q];
    
inverter_dxdt = [f;g];

end
