function f = inner_loop(vc_d_bar, vc_q_bar, vc, wc, ig_d, ig_q,  params)

  %get parameters 
  
  rv = params.rv;
  lv = params.lv;

  

 f = [
     
    % Eq. 21 d(xi_d)/dt=
    v_bar_d - eg_d;
    
    % Eq. 21 d(xi_q)/dt=
    v_bar_q - eg_q;
    
    % Eq. 23 d(gamma_d)/dt=
    is_bar_d-is_d;
    
    % Eq. 23 d(gamma_q)/dt=
    is_bar_q-is_q;
    
    % Alegbraic Eq's
    % Eq 20
    is_bar_d-(kp_v*(v_bar_d-eg_d)+...
    ki_v*xi_d - wc*cf*eg_q + kf_i*ig_d);
    
    is_bar_q-(kp_v*(v_bar_d-eg_q)+...
    ki_v*xi_q + wc*cf*eg_d + kf_i*ig_q);
    
    % Eq 21
    im_bar_d-(kp_i*(is_bar_d-is_d)+...
    ki_i*gamma_d - wc*lf*is_q + kf_v*eg_d);
    
    im_bar_q-(kp_i*(is_bar_q-is_q)+...
    ki_i*gamma_q + wc*lf*is_d + kf_v*eg_q);

    ];

end

