function ic = init_cond(p)
parameters;
vn_d = 1;
vn_q = 0;
ig_d = 0.1;
ig_q = 0;
epsilon=0;
eg_hat_d=1;
eg_hat_q=0;
pc_tilde=0.5;
qc_tilde=0;
theta_pll=2*pi*50;
theta_apc=2*pi*50;
pc=p;
qc=0;
vc=1.01;
wc = 1;
xi_d=0;
xi_q=0;
gamma_d=0;
gamma_q=0;


eg_d=(lf/wb)*((rf*wb/lf)*ig_d-wb*wc*ig_q)+vn_d;

eg_q=(lf/wb)*(wb*wc*ig_d+(rf*wb/lf)*ig_q)+vn_q;

is_d=(cf/wb)*(-wc*wb*eg_q)-ig_d;

is_q=(cf/wb)*(wc*wb*eg_d)-ig_q;

%vm_d=(lf/wb)*((rf*wb/lf)*is_d-wb*wc*is_q)+eg_d;

%vm_q=(lf/wb)*(wb*wc*is_d+(rf*wb/lf)*is_q)+eg_q;

vc_bar_d=vc-rv*ig-d+wc*lv*ig_q;

vc_bar_q=-rv*ig_q+wc*lv*ig_d;

is_bar_d=kp_v*(v_bar_d-eg_d)+ki_v*xi_d - wc*cf*eg_q + kf_i*ig_d;

is_bar_q=kp_v*(v_bar_d-eg_q)+ki_v*xi_q + wc*cf*eg_d + kf_i*ig_q;

vm_d=kp_i*(is_bar_d-is_bar_d)+ki_i*gamma_d - wc*lf*is_q + kf_v*eg_d;

vm_q=kp_i*(is_bar_q-is_bar_q)+ki_i*gamma_q + wc*lf*is_d + kf_v*eg_q;

ic = [  eg_d;
        eg_q;
        ig_d;
        ig_q;
        is_d;
        is_q;
        xi_d;
        xi_q;
        gamma_d;
        gamma_q;
        epsilon;
        theta_apc;
        theta_pll;
        pc_tilde;
        qc_tilde;
        pc;
        qc;
        eg_hat_d;
        eg_hat_q;
        wpll;
        wc;
        vc;
        vc_bar_d;
        vc_bar_q;
        is_bar_d;
        is_bar_q;
        im_bar_d;
        im_bar_q;
        vm_d;
        vm_q];
end
