function f = PLL(eg_d, eg_q, eg_hat, squig_pll, epsilon, w_pll, wc, params)

   %get parameters
   
    wb = params.wb; % PLL filter in rad/s
    kp_pll = params.kp_pll; % PLL proportional gain
    ki_pll = params.ki_pll; % PLL integraql gain 

f = [
      % Eq. 8 d(squig_pll)/dt = 
      (w_pll-wc)*wb;
      
      % Eq. 10 d(epsilon)/dt = 
      imag(eg_hat);
      
      % Eq. 11 d(theta_pll)/dt=
      w_pll*wb;
      
      % Algebraic Eq's
      % Eq. 7 eg_hat
      eg_hat-(eg_d+1i*eg_q)*exp(-1i*squig_pll);
      
      % Eq. 9 
      wpll-wo-kp_pll*imag(eg_hat)-ki_pll*epsilon
      
       ];
end

