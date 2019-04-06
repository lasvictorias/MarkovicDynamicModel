function f = outer_control_loop(wc, vc, pc_tilde, pc, qc_tilde, qc, params)

  %get parameters 
  
  fz = params.fz;
  pc_ref = params.pc_ref;
  rp = params.rp;
  

 f = [
     
    % Eq 13. d(pc_tilde)/dt = 
    2*pi*fz*(pc - pc_tilde);
    
    % Eq. 14 d(wc)/dt = 
    (1/(2*pi*fz*rp))*(pc_ref - pc)-...
    (1/(2*pi*fz*rp))*(1/rp)*(wc-wc-ref);
    
    % Eq. 17 d(qc_tilde)/dt=
    2*pi*fz*(qc-qc_tilde);

    % Alegbraic Eq's
    % Eq 12
    wc-(wc_ref+rp*(pc_ref - pc_est));
    
    % Eq 16 
    vc - (vc_ref + rq(qc_ref - qc-tilde));

    ];

end

