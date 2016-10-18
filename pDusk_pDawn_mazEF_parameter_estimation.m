%% pDusk_pDawn_mazEF_parameter_estimation
% Look for the following tags: %!remove/ %!missing/ %!adjust/ %!..        %!remove
% Description      %!adjust


%% Clear History
clc;                % clear display
clear all;          % clear Workspace
close all;          % close all figures 
clear workspace;    % clear workspace

tic;
%% 
%% Common Parameters of pDusk and pDawn
time_interval = 0.010;      % define time steps for the figures. 0.01h means taking a time point every 3.6s.
n_parameters = 10;          % n_parameters = number of additional parameters to implement mazEF
n_pDusk = 23;
n_pDawn = 28;
n = 10000;                  % number of the sample parameter sets (for reading the data)
System = 'pDusk_pDawn';     % for writing data files
Underscore = '_';           % for writing data files


%% Unit Conversion of Light Intensities
% As the ODE system follows the units µmol, meters  (m) and hours (h), the
% given light intensities from Ohlendorf et al. (2012) in [µW cm^-2] needed
% to be transformed. This was done according to: Berthold Technologies 
% (https://www.berthold.com/en/bio/how-do-i-convert-irradiance-photon-flux)
% Addition: There is a correction to the given Ohlendorf nW cm^-2 given. 
% It should read µW cm^-2 (Corrigendum: DOI: 10.1016/j.jmb.2013.11.005).

% Constants
h = 6.63*10^-34;                        % [Js] Planck's constant
c = 2.988*10^8;                         % [m/s] speed of light
lambda_nm = 470;                        % [nm] wavelength
N_A = 6.022*10^23;                      % [mol^-1] Avogadro's constant


%% Reading the Generated Data
Underscore = '_';

% Output_pDusk
spec = '%d%s%.3f%s%s%sOutput_pDusk.csv'; 
str = sprintf(spec, n, Underscore, time_interval, Underscore, System, Underscore);
Output_pDusk = csvread(str);

% Output_pDawn
spec = '%d%s%.3f%s%s%sOutput_pDawn.csv'; 
str = sprintf(spec, n, Underscore, time_interval, Underscore, System, Underscore);
Output_pDawn = csvread(str);

% BestOutput_pDusk
spec = '%d%s%.3f%s%s%sBestOutput_pDusk.csv'; 
str = sprintf(spec, n, Underscore, time_interval, Underscore, System, Underscore);
BestOutput_pDusk = csvread(str);

% BestOutput_pDawn
spec = '%d%s%.3f%s%s%sBestOutput_pDawn.csv'; 
str = sprintf(spec, n, Underscore, time_interval, Underscore, System, Underscore);
BestOutput_pDawn = csvread(str);

% Score
spec = '%d%s%.3f%s%s%sScore.csv'; 
str = sprintf(spec, n, Underscore, time_interval, Underscore, System, Underscore);
Score = csvread(str);

% Parameter_Set
spec = '%d%s%.3f%s%s%sParameter_Set.csv'; 
str = sprintf(spec, n, Underscore, time_interval, Underscore, System, Underscore);
Parameter_Set = csvread(str);


%% Redefining n for generating parameter samples for mazEF addition
n = 1000;


%% Generating mazEF Parameter Sets
% 10 new parameters for implementing mazEF
lhs_set = lhsdesign(n*2.1, n_parameters);                   % generating parameter estimates between 0 and 1 using Latin Hypercube Sampling 
lhs = [(1:(n*2.1))' (-log(lhs_set(:,1:n_parameters)))];     % lognormal transformation for parameter 1:10

% adding constraint that degradation rate of mazF_p must be higher than degradation rate of mazE_p
index = [];
for i = 1:(n*2.1)
    if lhs(i,9) > lhs(i,3);
        index = [index i];
    end
end
lhs = lhs(index,:);
lhs = lhs(1:n,:);


%% Sorting the newly generated parameters into the previously obtained parameters from pDusk/ pDawn

% Reading the Best Score
BestScore = Score(Score(:,2) == min(Score(:,2)),:);

% Selecting the Parameters of the Best Parameter Set
best_p = Parameter_Set(BestScore(1),:);


best_p_matrix = [];
for i = 1:21
    I(1:n) = best_p(i);
    best_p_matrix = [best_p_matrix; I];
end
best_p_matrix = best_p_matrix';

parameter_set_pDusk = [best_p_matrix(:,1:13), lhs];
% ------------------------- From best_p
% k_1           = p(1)      production rate of y_DD
% k_2           = p(2)      relaxation rate (tau) of YF1 'buffer'-system
% k_3           = p(3)      conversion cross-section (sigma) of light-intensity activated production rate
% beta_1        = p(4)      degradation rate of y_DD
% beta_2        = p(5)      degradation rate of y_DL/LD
% beta_3        = p(6)      degradation rate of y_LL
% k_4           = p(7)      production rate of j_i
% k_5           = p(8)      spontaneous de-phosphorylation rate
% beta_4        = p(9)      degradation rate of j_i
% k_6           = p(10)     production rate of j_a depending on the concentration of y_DD and j_i
% beta_5        = p(11)     degradation rate of j_a
% V_max         = p(12)     V_max of cT_m production based on j_a (same as V_max of RFP_m production based on j_a)
% K_m           = p(13)     K_m of cT_m production based on j_a (same as K_m of RFP_m production based on j_a)
% ------------------------- New Parameter Estimates pDusk mazEF
% beta_6        = p(14)     degradation rate of e_m                                                 = p(22) in pDawn_function_const_mazF
% k_7           = p(15)     production rate from e_m to e_p                                         = p(23) in pDawn_function_const_mazF
% beta_7        = p(16)     degradation rate of e_p                                                 = p(24) in pDawn_function_const_mazF
% k_10          = p(17)     dissociation rate of complex ef (lumped/ simplified)                    = p(25) in pDawn_function_const_mazF
% k_11          = p(18)     rate of ef-complex formation (lumped/ simplified)                       = p(26) in pDawn_function_const_mazF
% k_12          = p(19)     production rate of f_m based on constitutive promoter                   = p(27) in pDawn_function_const_mazF
% beta_10       = p(20)     degradation rate of f_m                                                 = p(14) in pDawn_function_const_mazF
% k_13          = p(21)     production rate from f_m to f_p                                         = p(15) in pDawn_function_const_mazF
% beta_11       = p(22)     degradation rate of f_p                                                 = p(16) in pDawn_function_const_mazF
% beta_12       = p(23)     degradation of complex ef                                               = p(28) in pDawn_function_const_mazF

parameter_set_pDawn = [best_p_matrix(:,1:13),   lhs(:,7:9), best_p_matrix(:,17:21), lhs(:,1:6), lhs(:,10)];
%                      p(1:13),                 p(14:16),   p(17:21),                p(22:27);   p(28)
%                      common in pDusk/pDawn    beta_10     
%                                               k_13        
%                                               beta_11     
% ------------------------- From best_p
% k_1           = p(1)      production rate of y_DD
% k_2           = p(2)      relaxation rate (tau) of YF1 'buffer'-system
% k_3           = p(3)      conversion cross-section (sigma) of light-intensity activated production rate
% beta_1        = p(4)      degradation rate of y_DD
% beta_2        = p(5)      degradation rate of y_DL/LD
% beta_3        = p(6)      degradation rate of y_LL
% k_4           = p(7)      production rate of j_i
% k_5           = p(8)      spontaneous de-phosphorylation rate
% beta_4        = p(9)      degradation rate of j_i
% k_6           = p(10)     production rate of j_a depending on the concentration of y_DD and j_i
% beta_5        = p(11)     degradation rate of j_a
% V_max         = p(12)     V_max of cT_m production based on j_a (same as V_max of RFP_m production based on j_a)
% K_m           = p(13)     K_m of cT_m production based on j_a (same as K_m of RFP_m production based on j_a)
% ...           ...         ...
% beta_8        = p(17)     degradation rate of lambda phage inhibitor mRNA (cI_m)
% k_8           = p(18)     production rate of cI_p depending on cI_m
% beta_9        = p(19)     degradation rate of cI_p
% k_9           = p(20)     maximal production rate of RFP_m (maximal transcription rate of the promoter)
% K_d           = p(21)     dissociation constant of cI_p at RFP_m promoter
% ------------------------- New Parameter Estimates pDawn mazEF
% beta_10       = p(14)     degradation rate of f_m                                                 = p(20) in pDusk_function_const_mazF
% k_13          = p(15)     production rate from f_m to f_p                                         = p(21) in pDusk_function_const_mazF
% beta_11       = p(16)     degradation rate of f_p                                                 = p(22) in pDusk_function_const_mazF
% ...           ...         ...
% beta_6        = p(22)     degradation rate of e_m                                                 = p(14) in pDusk_function_const_mazF
% k_7           = p(23)     production rate from e_m to e_p                                         = p(15) in pDusk_function_const_mazF
% beta_7        = p(24)     degradation rate of e_p                                                 = p(16) in pDusk_function_const_mazF
% k_10          = p(25)     dissociation rate of complex ef (lumped/ simplified)                    = p(17) in pDusk_function_const_mazF
% k_11          = p(26)     rate of ef-complex formation (lumped/ simplified)                       = p(18) in pDusk_function_const_mazF
% k_12          = p(27)     production rate of e_m based on constitutive promoter                   = p(19) in pDusk_function_const_mazF
% beta_12       = p(28)     degradation of complex ef                                               = p(23) in pDusk_function_const_mazF


%% 
%% pDusk
% Defining variables for pDusk
t   = 0:time_interval:17;       % Ohlendorf et al. (2013) induced the bacteria for 17h, i.e. 61200s. Thus the same time frame is used
x0  = [1 0 0 0 1 0 0 0 0 0];    % initial guesses for the variables of the ODE system. In the beginning state, the system is assumed to be in dark conditions, thus y_DD = 100% and j_a = 100% ("all" inactive FixJ (j_i) is phosphorylated)
p = zeros(n_pDusk);             % creating a zero-filled parameter vector


%% Loading the Data Points of Ohlendorf et al. (2012) for pDusk
% The data points of figure 2a in Ohlendorf et al. (2012) were retrieved
% through WebPlotDigitizer (http://arohatgi.info/WebPlotDigitizer) by (Ankit
% Rohatgi) to be used for the scoring function to evaluate the parameter
% sets generated by the latin hypercube sampling.
DataPoints_pDusk = csvread('DataPoints_pDusk.csv');     % reading the data points of figure 2a from Ohlendorf et al. (2012)
lightintensities_pDusk = [DataPoints_pDusk(:,1)'];      
RFP_datapoints_pDusk = [DataPoints_pDusk(:,2)'];

% Conversion of Light Intensities
N_pDusk = lightintensities_pDusk*10^-2;                 % µW*cm^-2 --- *10^-2 ---> W*m^-2
E_p = h*c/lambda_nm;                                    % distinct energy quanta of a photon
N_p = N_pDusk*lambda_nm*1E6/(1.988*1E-16);              % number of photons; lambda in [nm]
E_q = N_p/N_A;                                          % photonflux N = E_q in [Âµmol*m^-2*s^-1]
E_q_h = E_q*60*60;                                      % photonflux in [Âµmol*m^-2*h^-1]
N_pDusk = E_q_h;


%% Simulating Each parameter_set_pDusk(j) per Each Light Intensity N_pDusk(i)
X_Output = [];
P = [];
N_all = [];
RFP_all = [];
for j = 1:n;
    for i = 1:length(N_pDusk);
        [T, X] = ode45(@(t,x) pDusk_function_const_mazF(t,x,parameter_set_pDusk(j,:),N_pDusk(i)), t, x0);
        X_Output = [X_Output; X(end,:)];
    end 
    J(1:length(N_pDusk)) = j;
    P = [P; J'];
    N_all = [N_all; N_pDusk'];
    RFP_all = [RFP_all; RFP_datapoints_pDusk'];
end
Output_pDusk = [P, N_all, RFP_all, X_Output];     
% Column 1      = parameter_set_pDusk(j,:)
% Column 2      = N_pDusk
% Column 3      = RFP datapoints (pDusk)
% Column 4:13   = y1:y10
% Column 10     = y(7) = mazE
% Column 12     = y(9) = mazF


%% Recording the mazE/mazF (A/T) Ratio
P = [];
for i = 1:n;
   P = [P; repmat(i, length(N_pDusk), 1)];
end
pDusk_Ratio = [P, repmat((1:length(N_pDusk))',n,1), (Output_pDusk(:,10)./Output_pDusk(:,12)), Output_pDusk(:,10), Output_pDusk(:,12)];


%%
%% pDawn
% Defining variables for pDawn
t   = 0:time_interval:17;   % Ohlendorf et al. (2013) induced the bacteria for 17h, i.e. 61200s. Thus the same time frame is used
x0  = [1 0 0 0 1 0 1 0 0 0 0 0];  % initial guesses for the variables of the ODE system. In the beginning state, the system is assumed to be in dark conditions, thus y_DD = 100% and j_a = 100% ("all" inactive FixJ (j_i) is phosphorylated)
p = zeros(n_pDawn);    % creating a zero-filled parameter vector


%% Loading the Data Points of Ohlendorf et al. (2012) for pDawn
% The data points of figure 2a in Ohlendorf et al. (2012) were retrieved
% through WebPlotDigitizer (http://arohatgi.info/WebPlotDigitizer) by (Ankit
% Rohatgi) to be used for the scoring function to evaluate the parameter
% sets generated by the latin hypercube sampling.
DataPoints_pDawn = csvread('DataPoints_pDawn.csv');     % reading the data points of figure 2a from Ohlendorf et al. (2012)
lightintensities_pDawn = [DataPoints_pDawn(:,1)'];      
RFP_datapoints_pDawn = [DataPoints_pDawn(:,2)'];

% Conversion of Light Intensities
N_pDawn = lightintensities_pDawn*10^-2;                 % ÂµW*cm^-2 --- *10^-2 ---> W*m^-2
E_p = h*c/lambda_nm;                                    % distinct energy quanta of a photon
N_p = N_pDawn*lambda_nm*1E6/(1.988*1E-16);              % number of photons; lambda in [nm]
E_q = N_p/N_A;                                          % photonflux N = E_q in [Âµmol*m^-2*s^-1]
E_q_h = E_q*60*60;                                      % photonflux in [Âµmol*m^-2*h^-1]
N_pDawn = E_q_h;


%% Simulating Each parameter_set_pDawn(j) per Each Light Intensity N_pDawn(i)
X_Output = [];
P = [];
N_all = [];
RFP_all = [];
for j = 1:n
    for i = 1:length(N_pDawn)
        [T, X] = ode45(@(t,x) pDawn_function_const_mazE(t,x,parameter_set_pDawn(j,:),N_pDawn(i)), t, x0);
        X_Output = [X_Output; X(end,:)];
    end 
    J(1:length(N_pDawn)) = j;
    P = [P; J'];
    N_all = [N_all; N_pDawn']; 
end
Output_pDawn = [P, N_all, X_Output];     
% Column 1      = parameter_set_pDawn(j,:)
% Column 2      = N_pDawn
% Column 3:14   = y1:y12
% Column 13     = y(11) = mazE
% Column 11     = y(9) = mazF


%% Recording the mazE/mazF (A/T) Ratio
P = [];
for i = 1:n;
   P = [P; repmat(i, length(N_pDawn), 1)];
end
pDawn_Ratio = [P, repmat((1:length(N_pDawn))',n,1), (Output_pDawn(:,13)./Output_pDawn(:,11)), Output_pDawn(:,13), Output_pDawn(:,11)];


%%
%% Common Analysis
% %% Distribution of mazE/mazF (A/T) Ratio in pDusk and pDawn
% %% An Average Ratio Between Antitoxin and Toxin
% % Fasani and Savageau (2013)
% % Cases     Condition                               Ratios                              Comment
% % 1-16      µ           > lambda_A  > lambda_T      A/T = sigma                         implying the growth rate is relatively high, which is expected during normal operation
% % 17-32     lambda_T    > µ         > lambda_A      -                                   but ?A is actually greater than ?T in every well-studied toxin?antitoxin system, so those cases are ignored
% % 33-48     lambda_A    > µ         > lambda_T      A/T = sigma*µ_max/lambda_A          implying growth is slower, but dilution of the toxin is still significant, a situation we describe as transitional
% % 49-64     lambda_A    > lambda_T  > µ             A/T = sigma*lambda_T/ lambda_A      indicating the growth rate is near zero and described as static operation
% 
% % % From Fasani and Savageau (2013) Supplementary Material:
% % I assume that case 49-64 is present, thus: 
% % lambda_A    > lambda_T  > µ 
% % A/T = sigma*lambda_T/ lambda_A
% 
% % sigma     = 10                                    (assuming it is near the value they chose for their model)
% % lambda_T  = log(2)/ 4 [h]      = 0.1733 [1/h]
% % lambda_A  = log(2)/ 0.5 [h]    = 1.3863 [1/h]
% 
% AT = 10*(log(2)/4)/(log(2)/ 0.5);
% 
% 
% %% Constraining the A/T ratios
% constrained_pDusk_Ratio = pDusk_Ratio(pDusk_Ratio(:,3) > 1.1 & pDusk_Ratio(:,3) < 1.4 & pDusk_Ratio(:,2) == 1,:);
% constrained_pDawn_Ratio = pDawn_Ratio(pDawn_Ratio(:,3) > 10 & pDawn_Ratio(:,2) == 1,:);
% 
% combined_constrains = [];
% pDusk_ratio = [];
% for i = 1:size(constrained_pDusk_Ratio)
%     combined_constrains = [combined_constrains; constrained_pDawn_Ratio(constrained_pDawn_Ratio(:,1) == constrained_pDusk_Ratio(i,1), :)];
%     for k = 1:length(constrained_pDawn_Ratio(:,1));
%         if constrained_pDawn_Ratio(k,1) == constrained_pDusk_Ratio(i,1)
%             pDusk_ratio = [pDusk_ratio, constrained_pDusk_Ratio(i,3)];
%         end
%     end
% end
% combined_constrains = [combined_constrains pDusk_ratio'];
% 
% 
% %% Searching for the Best Parameter Set for pDusk and pDawn
% min_AT_difference = min((constrained_pDusk_Ratio(:,3)-AT).^2);
% if combined_constrains(:,4) == min_AT_difference;
%     best_combined_constrains = combined_constrains(combined_constrains(:,4) == min_AT_difference, :);
% else
%     min_pDawn_ratio = min(combined_constrains(:,3));
%     best_combined_constrains = combined_constrains(combined_constrains(:,3) == min_pDawn_ratio, :);
% end
% 
% 
% %% Storing the best parameter sets per pDusk/ pDawn
% best_p_pDusk = parameter_set_pDusk(best_combined_constrains(1),:);
% best_p_pDawn = parameter_set_pDawn(best_combined_constrains(1),:);


%% Write Outputs to csv files
% Parameter Set pDusk
spec = '%d%s%.3f%s%s%sParameter_Set_pDusk.csv'; 
str = sprintf(spec, n, Underscore, time_interval, Underscore, System, Underscore);
dlmwrite(str, parameter_set_pDusk);

% Parameter Set pDawn
spec = '%d%s%.3f%s%s%sParameter_Set_pDawn.csv'; 
str = sprintf(spec, n, Underscore, time_interval, Underscore, System, Underscore);
dlmwrite(str, parameter_set_pDawn);

% Output_pDusk_mazEF
spec = '%d%s%.3f%s%s%sOutput_pDusk_mazEF.csv'; 
str = sprintf(spec, n, Underscore, time_interval, Underscore, System, Underscore);
dlmwrite(str, Output_pDusk);

% Output_pDawn_mazEF
spec = '%d%s%.3f%s%s%sOutput_pDawn_mazEF.csv'; 
str = sprintf(spec, n, Underscore, time_interval, Underscore, System, Underscore);
dlmwrite(str, Output_pDawn);

% pDusk Ratio
spec = '%d%s%.3f%s%s%spDusk_Ratio.csv'; 
str = sprintf(spec, n, Underscore, time_interval, Underscore, System, Underscore);
dlmwrite(str, pDusk_Ratio);

% pDawn Ratio
spec = '%d%s%.3f%s%s%spDawn_Ratio.csv'; 
str = sprintf(spec, n, Underscore, time_interval, Underscore, System, Underscore);
dlmwrite(str, pDawn_Ratio);

toc;