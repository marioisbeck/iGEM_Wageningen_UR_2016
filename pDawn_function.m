function dydt = pDawn_function(t,y,p,N)
%% Ordinary Differential Equations Describing the pDusk System
dydt = zeros(size(y));      % creating the vector y 

%% Variables
% y_DD          = y(1)          YF1 homodimer in dark-dark state
% y_DL/LD       = y(2)          lumped YF1 homodimer in both dark-light (DL) and light-dark (LD) state
% y_LL          = y(3)          YF1 homodimer in light-light state
% j_i           = y(4)          inactive form of FixJ (mRNA stage of FixJ is lumped)
% j_a           = y(5)          active form of FixJ
% cI_m          = y(6)          lambda phage inhibitor mRNA 
% cI_p          = y(7)          lambda phage inhibitor protein 
% RFP_m         = y(6 or 8)     mRNA form of Red Fluorescent Protein (RFP)
% RFP_p         = y(7 or 9)     protein form of RFP


%% Parameters
% k_1           = p(1)      	production rate of y_DD
% k_2           = p(2)      	relaxation rate of YF1 based on tau
% k_3           = p(3)      	conversion cross-section (sigma) of light-intensity activated production rate
% beta_1        = p(4)      	degradation rate of y_DD
% beta_2        = p(5)      	degradation rate of y_DL/LD
% beta_3        = p(6)      	degradation rate of y_LL
% k_4           = p(7)      	production rate of j_i
% k_5           = p(8)      	spontaneous de-phosphorylation rate
% beta_4        = p(9)      	degradation rate of j_i
% k_6           = p(10)     	production rate of j_a depending on the concentration of y_DD and j_i
% beta_5        = p(11)    		degradation rate of j_a
% V_max         = p(12)     	V_max of cI_m production based on j_a (same as V_max of RFP_m production based on j_a)
% K_m           = p(13)     	K_m of cI_m production based on j_a (same as K_m of RFP_m production based on j_a)
% beta_6        = p(14)     	degradation rate of RFP_m
% k_7           = p(15)     	translation rate from RFP_m to RFP_p
% beta_7        = p(16)     	degradation rate of RFP_p

%% Unique Parameters in pDawn
% beta_8        = p(17)     	degradation rate of lambda phage inhibitor mRNA
% k_8           = p(18)     	production rate of cI_p depending on cI_m
% beta_9        = p(19)     	degradation rate of cI_p
% k_9           = p(20)     	maximal production rate of RFP_m (maximal transcription rate of the promoter)
% K_d           = p(21)     	dissociation constant of cI_p at RFP_m promoter


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

% dcI_m/dt      = ((V_max*j_a) /(K_m+j_a))      - beta_8*cI_m;
dydt(6)         = ((p(12)*y(5))/(p(13)+y(5)))   - p(17)*y(6);

% dcI_p/dt      = k_8  *cI_m - beta_9*cI_p;
dydt(7)         = p(18)*y(6) - p(19)*y(7);

% dRFP_m/dt     = k_9  *(1/(1 + (cI_p/K_d)^2))   - beta_6*RFP_m;
dydt(8)         = p(20)*(1/(1 + (y(7)/p(21))^2)) - p(14)*y(8);  % http://parts.igem.org/wiki/index.php?title=Part:BBa_R0051 -> two binding sites for cI -> Hill Coefficient 2

% dRFP_p/dt     = k_7*RFP_m - beta_7*RFP_p;
dydt(9)         = p(15)*y(8) - p(16)*y(9);


end

