function dydt = pDawn_function_const_mazE(t,y,p,N)
%% Ordinary Differential Equations Describing the pDusk System
dydt = zeros(size(y));      % creating the vector y


%% Variables
% y_DD          = y(1)      YF1 homodimer in dark-dark state (y_DD)
% y_DL/LD       = y(2)      lumped YF1 homodimer in both dark-light (DL) and light-dark (LD) state (y_DL/LD)
% y_LL          = y(3)      YF1 homodimer in light-light (y_LL) state
% j_i           = y(4)      inactive form of FixJ (j_i) [mRNA stage of FixJ is lumped]
% j_a           = y(5)      active form of FixJ (j_a)
% cI_m          = y(6)      lambda phage inhibitor mRNA (cI_m)
% cI_p          = y(7)      lambda phage inhibitor protein (cI_p)
% f_m           = y(8)     	mRNA form of mazF (f_m)
% f_p           = y(9) 		protein form of mazF (f_p)
% e_m 			= y(10) 	mRNA form of mazE (e_m)
% e_p 			= y(11) 	protein form of mazE (e_p)
% ef 			= y(12) 	inactive complex form of mazE-mazF (ef) [lumped/ simplified complex formation]


%% Parameters
% k_1           = p(1)      production rate of y_DD
% k_2           = p(2)      relaxation rate (tau) of YF1 'buffer'-system        %!adjust
% k_3           = p(3)      conversion cross-section (sigma) of light-intensity activated production rate
% beta_1        = p(4)      degradation rate of y_DD
% beta_2        = p(5)      degradation rate of y_DL/LD
% beta_3        = p(6)      degradation rate of y_LL
% k_4           = p(7)      production rate of j_i
% k_5           = p(8)      spontaneous de-phosphorylation rate
% beta_4        = p(9)      degradation rate of j_i
% k_6           = p(10)     production rate of j_a depending on the concentration of y_DD and j_i
% beta_5        = p(11)     degradation rate of j_a
% V_max         = p(12)     V_max of cI_m production based on j_a (same as V_max of f_m production based on j_a)
% K_m           = p(13)     K_m of cI_m production based on j_a (same as K_m of f_m production based on j_a)

%% Varies to pDusk_function_const_mazF
% beta_12       = p(14)     degradation rate of f_m 												= p(20) in pDusk_function_const_mazF
% k_13          = p(15)     production rate from f_m to f_p 										= p(21) in pDusk_function_const_mazF
% beta_13       = p(16)     degradation rate of f_p 												= p(22) in pDusk_function_const_mazF
% ... 			... 		...
% beta_10       = p(22)     degradation rate of e_m 												= p(14) in pDusk_function_const_mazF
% k_14          = p(23)     production rate from e_m to e_p 										= p(15) in pDusk_function_const_mazF
% beta_11       = p(24)     degradation rate of e_p 												= p(16) in pDusk_function_const_mazF

%% Added Parameters mazEF
% k_10			= p(25)  	dissociation rate of complex ef (lumped/ simplified) 					= p(17) in pDusk_function_const_mazF
% k_11 			= p(26)  	rate of ef-complex formation (lumped/ simplified) 						= p(18) in pDusk_function_const_mazF
% k_12 			= p(27)  	production rate of e_m based on constitutive promoter 				 	= p(19) in pDusk_function_const_mazF
% beta_14 		= p(28)  	degradation of complex ef 												= p(23) in pDusk_function_const_mazF

%% Unique Parameters in pDawn
% beta_8        = p(17)     degradation rate of lambda phage inhibitor mRNA (cI_m)
% k_8           = p(18)     production rate of cI_p depending on cI_m
% beta_9        = p(19)     degradation rate of cI_p
% k_9           = p(20)     maximal production rate of f_m (maximal production rate of the promoter)
% K_d           = p(21)     dissociation constant of cI_p at f_m promoter



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

% dcI_m/dt      = ((V_max*j_a) /(K_m+j_a))      - beta_8*cI_m
dydt(6)         = ((p(12)*y(5))/(p(13)+y(5)))   - p(17)*y(6);

% dcI_p/dt      = k_8  *cI_m - beta_9*cI_p
dydt(7)         = p(18)*y(6) - p(19)*y(7);

% df_m/dt       = k_9  *(1/(1 + (cI_p/K_d)^2))   - beta_12*f_m
dydt(8)         = p(20)*(1/(1 + (y(7)/p(21))^2)) - p(14)*y(8);  % http://parts.igem.org/wiki/index.php?title=Part:BBa_R0051 -> two binding sites for cI -> Hill Coefficient 2

% df_p/dt       = k_13*f_m    + k_10 *(4*ef)    - k_11 *(e_p^2)  *(4*(f_p^4))  - beta_13*f_p;
dydt(9)         = p(15)*y(8)  + p(17)*(4*y(12)) - p(18)*(y(11)^2)*(4*(y(9)^4)) - p(16)*y(9); % adjustment 2* added in front of y(10) and y(7) as well as ^2 and ^4 [according to Sontag E.D. (2015) - Mathematical Systems Biology, p.85]

% de_m/dt 		= k_12  - beta_10*e_m;
dydt(10) 		= p(27) - p(22)*y(10);

% de_p/dt 		= k_14*e_m     + k_10 *(2*ef)   	- k_11 *(2*(e_p^2))  *(f_p^4)  - beta_11*e_p;
dydt(11) 		= p(23)*y(10) + p(17)*(2*y(12)) - p(18)*(2*(y(11)^2))*(y(9)^4) - p(24)*y(11); % adjustment 4* added in front of y(10) and y(9) as well as ^2 and ^4 [according to Sontag E.D. (2015) - Mathematical Systems Biology, p.85]

% def/dt 		= k_11 *(2*e_p)  *(4*f_p)  - k_10*ef     - beta_14*ef;
dydt(12) 		= p(18)*(y(11)^2)*(y(9)^4) - p(17)*y(12) - p(28)*y(12);


end

