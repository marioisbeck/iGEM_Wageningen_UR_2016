%% pDusk_function
function dydt = pDusk_function_const_mazF(t,y,p,N)
%% Ordinary Differential Equations Describing the pDusk System
dydt = zeros(size(y));      % creating the vector y


%% Variables
% y_DD          = y(1)      YF1 homodimer in dark-dark state (y_DD)
% y_DL/LD       = y(2)      lumped YF1 homodimer in both dark-light (DL) and light-dark (LD) state (y_DL/LD)
% y_LL          = y(3)      YF1 homodimer in light-light (y_LL) state
% j_i           = y(4)      inactive form of FixJ (j_i) [mRNA stage of FixJ is lumped]
% j_a           = y(5)      active form of FixJ (j_a)
% e_m         	= y(6)     	mRNA form of mazE (e_m)
% e_p         	= y(7)		protein form of mazE (e_p)
% f_m 			= y(8) 		mRNA form of mazF (f_m)
% f_p 			= y(9)		protein form of mazF (f_p)
% ef 			= y(10) 	inactive complex form of mazE-mazF (ef) [lumped/ simplified complex formation]


%% Parameters
% k_1           = p(1)      production rate of y_DD
% k_2           = p(2)      relaxation rate (tau) of YF1 
% k_3           = p(3)      conversion cross-section (sigma) of light-intensity activated production rate
% beta_1        = p(4)      degradation rate of y_DD
% beta_2        = p(5)      degradation rate of y_DL/LD
% beta_3        = p(6)      degradation rate of y_LL
% k_4           = p(7)      production rate of j_i
% k_5           = p(8)      spontaneous de-phosphorylation rate
% beta_4        = p(9)      degradation rate of j_i
% k_6           = p(10)     production rate of j_a depending on the concentration of y_DD and j_i
% beta_5        = p(11)     degradation rate of j_a
% V_max         = p(12)     V_max of e_m production based on j_a
% K_m           = p(13)     K_m of e_m production based on j_a

%% Varies to pDawn_function_const_mazF
% beta_10       = p(14)     degradation rate of e_m 			 									= p(22) in pDawn_function_const_mazF
% k_14          = p(15)     production rate from e_m to e_p 	 									= p(23) in pDawn_function_const_mazF
% beta_11       = p(16)     degradation rate of e_p 												= p(24) in pDawn_function_const_mazF

%% Added Parameters mazEF
% k_10			= p(17) 	dissociation rate of complex ef (lumped/ simplified) 					= p(25) in pDawn_function_const_mazF
% k_11 			= p(18) 	rate of ef-complex formation (lumped/ simplified) 						= p(26) in pDawn_function_const_mazF
% k_12 			= p(19) 	production rate of f_m based on constitutive promoter  				 	= p(27) in pDawn_function_const_mazF
% beta_12		= p(20) 	degradation rate of f_m 												= p(14) in pDawn_function_const_mazF
% k_13 			= p(21) 	production rate from f_m to f_p 				 						= p(15) in pDawn_function_const_mazF
% beta_13 		= p(22)		degradation rate of f_p 												= p(16) in pDawn_function_const_mazF
% beta_14 		= p(23) 	degradation of complex ef 												= p(28) in pDawn_function_const_mazF



%% Ordinary Differential Equations
% dy_DD/dt      = k_1  + 2*k_2*y_DL/LD  - 2*(N*k_3)*y_DD  - beta_1*y_DD;
dydt(1)         = p(1) + 2*p(2)*y(2)    - 2*(N*p(3))*y(1) - p(4)*y(1);

% dy_DL_LD/dt   = 2*(N*k_3)*y_DD  + 2*k_2*y_LL  - 2*k_2*y_DL/LD - 2*(N*k_3)*y_DL/LD  - beta_2*y_DL/LD;
dydt(2)         = 2*(N*p(3))*y(1) + 2*p(2)*y(3) - 2*p(2)*y(2)   - 2*(N*p(3))*y(2)    - p(5)*y(2);

% dy_LL/dt      = 2*(N*k_3)*y_DL/LD - 2*k_2*y_LL  - beta_3*y_LL;                    
dydt(3)         = 2*(N*p(3))*y(2)   - 2*p(2)*y(3) - p(6)*y(3);

% dj_i/dt       = k_4  + k_5*j_a   - beta_4*j_i;
dydt(4)         = p(7) + p(8)*y(5) - p(9)*y(4);

% dj_a/dt       = k_6*y_DD*j_i    -  beta_5*j_a;
dydt(5)         = p(10)*y(1)*y(4) - p(11)*y(5);

% de_m/dt       = ((V_max_e_p*j_a)/(K_m_e_p+j_a)) - beta_10*e_m;
dydt(6)         = ((p(12)*y(5))   /(p(13)+y(5)))  - p(14)*y(6);

% de_p/dt       = k_14*e_m    + k_10 *(2*ef)    - k_11 *(2*(e_p^2)) *(f_p^4)  - beta_11*e_p;
dydt(7)         = p(15)*y(6) + p(17)*(2*y(10)) - p(18)*(2*(y(7)^2))*(y(9)^4) - p(16)*y(7); % adjustment 2* added in front of y(10) and y(7) as well as ^2 and ^4 [according to Sontag E.D. (2015) - Mathematical Systems Biology, p.85]

% df_m/dt 		= k_12  - beta_12*f_m;
dydt(8) 		= p(19) - p(20)*y(8);

% df_p/dt 		= k_13*f_m   + k_10 *(4*ef)    - k_11 *(e_p^2)    *(4*(f_p^4))  - beta_13*f_p;
dydt(9) 		= p(21)*y(8) + p(17)*(4*y(10)) - p(18)*(y(7)^2)   *(4*(y(9)^4)) - p(22)*y(9); % adjustment 4* added in front of y(10) and y(9) as well as ^2 and ^4 [according to Sontag E.D. (2015) - Mathematical Systems Biology, p.85]

% def/dt 		= k_11 *(e_p^2) *(f_p^4)  - k_10*ef     - beta_14*ef;
dydt(10) 		= p(18)*(y(7)^2)*(y(9)^4) - p(17)*y(10) - p(23)*y(10);


end

