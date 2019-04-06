function f = virtual_impedance(wc, vc, pc_tilde, pc, qc_tilde, qc, params)

  %get parameters 
  
  rv = params.rv;
  pc_ref = params.pc_ref;
  rp = params.rp;
  

 f = [
    % Alegbraic Eq's
    % Eq 18
    vc_d_bar -(vc-rv*ig-d+wc*lv*ig_q);
    
    % Eq 19
    vc_d_q - (-rv*ig_q+wc*lv*ig_d);

    ];

end

