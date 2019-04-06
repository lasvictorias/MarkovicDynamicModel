function f = elec_interface(is_d, is_q, vm_d, vm_q, eg_d, eg_q, ig_d, ig_q, pc, qc, params)

 %get parameters
 
    wb = params.wb; 
    wc = params.wc; 
    lf = params.lf; % Voltage controller gain
    rf = params.rf; 
    cf = params.cf; 
    vn_d = params.vn_d;
    vn_q = params.vn_q;
    
f = [
    
     %d(is_d)/dt = 
      (wb/lf)*(vm_d-eg_d)-...
      (rf*wb/lf)*is_d+...
      wb*wc*is_q;
      
     %d(is_q)/dt = 
      (wb/lf)*(vm_q-eg_q)-...
      wb*wc*is_d-...
      (rf*wb/lf)*is_q;

      %d(ig_d)/dt = 
      (wb/lf)*(eg_d-vn_d)-...
      (rf*wb/lf)*ig_d+...
      wb*wc*ig_q;
      
     %d(ig_q)/dt = 
      (wb/lf)*(eg_q-vn_q)-...
      wb*wc*ig_d-...
      (rf*wb/lf)*ig_q;
      
      %d(eg_d)/dt =
      (wb/cf)*(is_d-ig_d)+...
      wc*wb*eg-q;
      
      %d(eg_q)/dt = 
      (wb/cf)*(is_q-ig_q)-...
      wc*wb*eg_d;
      
      %Alegrraic Eq's
      %P calc
      pc-...
      (eg_d*ig_d+...
      eg_q*ig_q);
      
      %Q calc
      qc-...
      (eg_q*ig_d-...
      ed_d*ig_q);

      ];
      
end

