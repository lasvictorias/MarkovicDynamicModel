    %% Differential equations
    % Low-pass filter 
    f = [f; p.vscF.wb/p.vscF.cf*y(Ns+s.vscF.istd) - p.vscF.wb/p.vscF.cf*y(Ns+s.vscF.igtd) + p.net.w*p.vscF.wb*y(Ns+s.vscF.egtq)];
    f = [f; p.vscF.wb/p.vscF.cf*y(Ns+s.vscF.istq) - p.vscF.wb/p.vscF.cf*y(Ns+s.vscF.igtq) - p.net.w*p.vscF.wb*y(Ns+s.vscF.egtd)];

    f = [f; p.vscF.wb/p.vscF.lf*y(Na+s.vscF.vmd) - p.vscF.wb/p.vscF.lf*y(Ns+s.vscF.egtd) - p.vscF.rf*p.vscF.wb/p.vscF.lf*y(Ns+s.vscF.istd) + p.vscF.wb*p.net.w*y(Ns+s.vscF.istq)]; 
    f = [f; p.vscF.wb/p.vscF.lf*y(Na+s.vscF.vmq) - p.vscF.wb/p.vscF.lf*y(Ns+s.vscF.egtq) - p.vscF.rf*p.vscF.wb/p.vscF.lf*y(Ns+s.vscF.istq) - p.vscF.wb*p.net.w*y(Ns+s.vscF.istd)]; 
    
    % Transformer
    f = [f; p.vscF.wb/p.vscF.lt*y(Ns+s.vscF.egtd) - p.vscF.wb/p.vscF.lt*y(Na+s.vscF.vtd) - p.vscF.rt*p.vscF.wb/p.vscF.lt*y(Ns+s.vscF.igtd) + p.vscF.wb*p.net.w*y(Ns+s.vscF.igtq)]; 
    f = [f; p.vscF.wb/p.vscF.lt*y(Ns+s.vscF.egtq) - p.vscF.wb/p.vscF.lt*y(Na+s.vscF.vtq) - p.vscF.rt*p.vscF.wb/p.vscF.lt*y(Ns+s.vscF.igtq) - p.vscF.wb*p.net.w*y(Ns+s.vscF.igtd)]; 
    
    % Active power control
    f =  [f; p.vscF.wf*(y(Na+s.vscF.pm) - y(Ns+s.vscF.pf))]; 
    f =  [f; y(Na+s.vscF.dwapc)*p.vscF.wb]; 

    % Reactive power control 
    f = [f; p.vscF.wf*(y(Na+s.vscF.qm)-y(Ns+s.vscF.qf))]; 

    % SRF voltage control        
    f = [f; y(Na+s.vscF.v0d_ref) - y(Na+s.vscF.egd)]; 
    f = [f; y(Na+s.vscF.v0q_ref) - y(Na+s.vscF.egq)];

    % SRF current control 
    f = [f; y(Na+s.vscF.isd_ref) - y(Na+s.vscF.isd)]; 
    f = [f; y(Na+s.vscF.isq_ref) - y(Na+s.vscF.isq)]; 
    
    %% Algebraic equations
    % dq-transformation (reference SRF change: network -> vsc)!
    g = [g; - y(Na+s.vscF.egd) + y(Ns+s.vscF.egtd) * cos(-y(Ns+s.vscF.dOapc)) - y(Ns+s.vscF.egtq) * sin(-y(Ns+s.vscF.dOapc))];
    g = [g; - y(Na+s.vscF.egq) + y(Ns+s.vscF.egtd) * sin(-y(Ns+s.vscF.dOapc)) + y(Ns+s.vscF.egtq) * cos(-y(Ns+s.vscF.dOapc))];

    g = [g; - y(Na+s.vscF.igd) + y(Ns+s.vscF.igtd)* cos(-y(Ns+s.vscF.dOapc)) - y(Ns+s.vscF.igtq) * sin(-y(Ns+s.vscF.dOapc))];
    g = [g; - y(Na+s.vscF.igq) + y(Ns+s.vscF.igtd)* sin(-y(Ns+s.vscF.dOapc)) + y(Ns+s.vscF.igtq) * cos(-y(Ns+s.vscF.dOapc))];

    g = [g; - y(Na+s.vscF.isd) + y(Ns+s.vscF.istd)* cos(-y(Ns+s.vscF.dOapc)) - y(Ns+s.vscF.istq) * sin(-y(Ns+s.vscF.dOapc))];
    g = [g; - y(Na+s.vscF.isq) + y(Ns+s.vscF.istd)* sin(-y(Ns+s.vscF.dOapc)) + y(Ns+s.vscF.istq) * cos(-y(Ns+s.vscF.dOapc))];
        
    % Power calculation
    g = [g; - y(Na+s.vscF.pm) + y(Na+s.vscF.egd) * y(Na+s.vscF.igd) + y(Na+s.vscF.egq) * y(Na+s.vscF.igq)]; 
    g = [g; - y(Na+s.vscF.qm) - y(Na+s.vscF.egd) * y(Na+s.vscF.igq) + y(Na+s.vscF.egq) * y(Na+s.vscF.igd)]; 
    
    % Active power control
    g = [g; - y(Na+s.vscF.wapc) + p.vscF.wstar(i) + p.vscF.Dp*(p.vscF.pstar(i)-y(Ns+s.vscF.pf))];
    g = [g; - y(Na+s.vscF.dwapc) + y(Na+s.vscF.wapc) - p.net.w]; 

    % Reactive power control
    g = [g; - y(Na+s.vscF.v) + p.vscF.vstar(i) + p.vscF.Dq*(p.vscF.qstar(i)-y(Ns+s.vscF.qf))];

    % Virtual impedance
    g = [g; - y(Na+s.vscF.v0d_ref) + y(Na+s.vscF.v) - p.vscF.rv*y(Na+s.vscF.igd) + y(Na+s.vscF.wapc)*p.vscF.lv*y(Na+s.vscF.igq)];
    g = [g; - y(Na+s.vscF.v0q_ref) + 0*y(Na+s.vscF.v) - p.vscF.rv*y(Na+s.vscF.igq) - y(Na+s.vscF.wapc)*p.vscF.lv*y(Na+s.vscF.igd)];

    % SRF voltage control
    g = [g; - y(Na+s.vscF.isd_ref) + p.vscF.Kpv*(y(Na+s.vscF.v0d_ref) - y(Na+s.vscF.egd)) + p.vscF.Kiv*y(Ns+s.vscF.xid) - y(Na+s.vscF.wapc)*p.vscF.c1*y(Na+s.vscF.egq) + p.vscF.Kffi*y(Na+s.vscF.igd)]; 
    g = [g; - y(Na+s.vscF.isq_ref) + p.vscF.Kpv*(y(Na+s.vscF.v0q_ref) - y(Na+s.vscF.egq)) + p.vscF.Kiv*y(Ns+s.vscF.xiq) + y(Na+s.vscF.wapc)*p.vscF.c1*y(Na+s.vscF.egd) + p.vscF.Kffi*y(Na+s.vscF.igq)];

    % SRF current control
    g = [g; - y(Na+s.vscF.vmd_ref) + p.vscF.Kpc*(y(Na+s.vscF.isd_ref) - y(Na+s.vscF.isd))+ p.vscF.Kic*y(Ns+s.vscF.gamd) - y(Na+s.vscF.wapc)*p.vscF.l1*y(Na+s.vscF.isq) + p.vscF.Kffv*y(Na+s.vscF.egd)]; 
    g = [g; - y(Na+s.vscF.vmq_ref) + p.vscF.Kpc*(y(Na+s.vscF.isq_ref) - y(Na+s.vscF.isq))+ p.vscF.Kic*y(Ns+s.vscF.gamq) + y(Na+s.vscF.wapc)*p.vscF.l1*y(Na+s.vscF.isd) + p.vscF.Kffv*y(Na+s.vscF.egq)]; 

    % Switching unit
    g = [g; - y(Na+s.vscF.vmd) + y(Na+s.vscF.vmd_ref) * cos(y(Ns+s.vscF.dOapc)) - y(Na+s.vscF.vmq_ref) * sin(y(Ns+s.vscF.dOapc))]; 
    g = [g; - y(Na+s.vscF.vmq) + y(Na+s.vscF.vmd_ref) * sin(y(Ns+s.vscF.dOapc)) + y(Na+s.vscF.vmq_ref) * cos(y(Ns+s.vscF.dOapc))]; 