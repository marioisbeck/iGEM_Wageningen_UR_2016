%% pDusk_pDawn_parameter_estimation_weighted_sum
% Look for the following tags: %!remove/ %!missing/ %!adjust/ %!..        %!remove
% Description      %!adjust


%% Clear History
clc;                % clear display
clear all;          % clear Workspace
close all;          % close all figures 
clear workspace;    % clear workspace


%% 
%% Common Parameters of pDusk and pDawn
time_interval = 0.010;      % define time steps for the figures. 0.01h means taking a time point every 3.6s.
n_pDusk  = 16;              % number of parameters of the ODE system of pDusk
n_pDawn  = 21;              % number of parameters of the ODE system of pDawn
n_parameters = 21;          % n_parameters = n_pDawn (for better readability of code)
n = 10000;                  % number of the sample parameter sets to be generated
System = 'pDusk_pDawn';     % for writing data files
Underscore = '_';           % for writing data files


%% Unit Conversion of Light Intensities
% As the ODE system follows the units �mol, meters  (m) and hours (h), the
% given light intensities from Ohlendorf et al. (2012) in [�W cm^-2] needed
% to be transformed. This was done according to: Berthold Technologies 
% (https://www.berthold.com/en/bio/how-do-i-convert-irradiance-photon-flux)
% Addition: There is a correction to the given Ohlendorf nW cm^-2 given. 
% It should read �W cm^-2 (Corrigendum: DOI: 10.1016/j.jmb.2013.11.005).

% Constants
h = 6.63*10^-34;                        % [Js] Planck's constant
c = 2.988*10^8;                         % [m/s] speed of light
lambda_nm = 470;                        % [nm] wavelength
N_A = 6.022*10^23;                      % [mol^-1] Avogadro's constant


%%
%% Generating Parameter Sets

lhs_set = lhsdesign(n, n_parameters);                   % generating parameter estimates between 0 and 1 using Latin Hypercube Sampling 


%% Adjusting the Distribution of the Generated Parameters
% As the Latin Hypercube Sampling generates uniformly distributed
% parameter sets, the generated data needs to be transformed. 

% For those parameters, where there is a clear estimate � standard 
% diviation given, a normal distribution around that estimate is assumed.
% For parameters, where there is no (good) estimate known or the statistics
% are missing, a negative lognormal transformation is applied as this
% reflects the distribution of parameters in nature.

% Parameter Distributions
x = lhs_set(:,2);                                   % Relaxation rate (k_2)
tau = 5900 + sqrt(2)*25*erfinv(2*x-1);              % 5900 � 25 from M�glich et al. (2009)
y = lhs_set(:,3);                                   % Conversion Cross-Section (k_3)
sigma = 1000 + sqrt(2)*250*erfinv(2*y-1);           % estimated based on Rausenberger et al. (2009) and Klose et al. (2015)

param_a = (-log(lhs_set(:,1)));                     % lognormal transformation for parameter 1
param_b = sigma.*10^-6;                             % conversion from [%!adjust] --> []         
param_c = (log(2)*((tau.*(3600^-1)).^-1));          % transformation of tau to relaxation rate             
param_d = (-log(lhs_set(:,4:n_pDusk)));             % lognormal transformation for parameters 4:16
param_e = (-log(lhs_set(:,(n_pDusk+1):n_pDawn)));   % lognormal transformation for parameters 17:22
parameter_set       = [param_a, param_b, param_c, param_d, param_e];     
parameter_set_pDusk = [param_a, param_b, param_c, param_d];
parameter_set_pDawn = [param_a, param_b, param_c, param_d, param_e];


%% 
%% pDusk
% Defining variables for pDusk
t   = 0:time_interval:17;   % Ohlendorf et al. (2013) induced the bacteria for 17h, i.e. 61200s. Thus the same time frame is used
x0  = [1 0 0 0 1 0 0];      % initial guesses for the variables of the ODE system. In the beginning state, the system is assumed to be in dark conditions, thus y_DD = 100% and j_a = 100% ("all" inactive FixJ (j_i) is phosphorylated)
p = zeros(n_pDusk);         % creating a zero-filled parameter vector


%% Loading the Data Points of Ohlendorf et al. (2012) for pDusk
% The data points of figure 2a in Ohlendorf et al. (2012) were retrieved
% through WebPlotDigitizer (http://arohatgi.info/WebPlotDigitizer) by (Ankit
% Rohatgi) to be used for the scoring function to evaluate the parameter
% sets generated by the latin hypercube sampling.
DataPoints_pDusk = csvread('DataPoints_pDusk.csv');     % reading the data points of figure 2a from Ohlendorf et al. (2012)
lightintensities_pDusk = [DataPoints_pDusk(:,1)'];      
RFP_datapoints_pDusk = [DataPoints_pDusk(:,2)'];

% Conversion of Light Intensities
N_pDusk = lightintensities_pDusk*10^-2;                 % �W*cm^-2 --- *10^-2 ---> W*m^-2
E_p = h*c/lambda_nm;                                    % distinct energy quanta of a photon
N_p = N_pDusk*lambda_nm*1E6/(1.988*1E-16);              % number of photons; lambda in [nm]
E_q = N_p/N_A;                                          % photonflux N = E_q in [�mol*m^-2*s^-1]
E_q_h = E_q*60*60;                                      % photonflux in [�mol*m^-2*h^-1]
N_pDusk = E_q_h;


%% Simulating Each parameter_set_pDusk(j) per Each Light Intensity N_pDusk(i)
X_Output = [];
P = [];
N_all = [];
RFP_all = [];
for j = 1:n;
    for i = 1:length(N_pDusk);
        [T, X] = ode45(@(t,x) pDusk_function(t,x,parameter_set_pDusk(j,:),N_pDusk(i)), t, x0);
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
% Column 4:10   = y1:y7
% Column 11     = normalized RFP_p (pDusk)
% Column 12     = StDev_RFP (pDusk)


%% Normalizing the Generated RFP_p Output (pDusk)
% Here the generated RFP_p concentration is normlized according to the highest light intensity readout of RFP in Ohlendorf et al. (2012)
% level (FL/ OD600) given in Ohlendorf et al. (2012) figure 2a.
est_RFP_level = [];         % estimated RFP level
K = Output_pDusk(:,1)*length(N_pDusk); 
for h = 1:length(Output_pDusk(:,1));
    est_RFP_level = [est_RFP_level; (0.032614/Output_pDusk(K(h),10))*Output_pDusk(h,10)]; 
    % 10                          = y(7)      = RFP_p
    % 0.032614                    = RFP datapoint from Ohlendorf et al. (2012) at N~99 (highest light intensity point) �W*cm^-2
    % Output_pDusk(K(h),10))      = taking the computed value of RFP_p = y(7) = Output(:,10) at the highest light intensity ~N99 of each parameter set
    % Output_pDusk(h,10)          = taking the computed value from each light intensity of each parameter set
end            
Output_pDusk = [Output_pDusk, est_RFP_level];
% Column 11     = normalized RFP_p (pDusk) 


%% Loading the Standard Deviation Values from Ohlendorf et al. (2012)
F = csvread('Standard_Deviation_Ohlendorf_pDusk.csv');
StDev_Ohlendorf_pDusk = [F(:,2)'];                                                            % this is StDev in logscale
StDev_Sorted_pDusk = StDev_Ohlendorf_pDusk([1,3,5,7,9,11,13,15]);                             % subsetting upper bound StDev values
StDev_Sorted_pDusk = [StDev_Sorted_pDusk; StDev_Ohlendorf_pDusk([2,4,6,8,10,12,14,16])];      % subsetting lower bound StDev values
% StDev_RFP_pDusk = (StDev_Sorted_pDusk(1,:)-StDev_Sorted_pDusk(2,:))/2;                      % extracting StDev values (not considering log scale)
StDev_RFP_pDusk = (exp(StDev_Sorted_pDusk(1,:))-exp(StDev_Sorted_pDusk(2,:)))/2;              % extracting StDev values (considering log scale)
%!question: I think I need to transform the StDev from a log scale to a normal scale?!

% Adding the StDev_RFP to the matrix "Output"
Output_pDusk = [Output_pDusk, repmat(StDev_RFP_pDusk,1,j)'];
% Column 12     = StDev_RFP_pDusk


%% Evaluating the Parameter Sets - The Scoring Function (pDusk)
% In order to understand whether a parameter set describes the given data
% by Ohlendorf et al. (2012) can be done by computing a scoring function.
% The idea is similar to the sum of squared errors (SSEs) approach and can
% be followed in more detail in Raue et al. (2009) equation IV.

% The Scoring Function
% The agreement of experimental data with the observables predicted by the model is measured by an objective function, commonly the weighted sum of squared residuals - Raue et al. (2009)
psetlength = [1:j]*length(N_pDusk);
H = [];
S = [];
for j = 1:n;
    if j == 1;
        for i = 1:psetlength(1);
            H = [H; ((Output_pDusk(i,3)-Output_pDusk(i,end-1))/StDev_RFP_pDusk(i))^2];
            % time point is               = 17h (as always the end-point readout at 17h of the computation is store in Output
            % Output_pDusk(i,3)           = datapoints (RFP_p readout [FL/ OD600]) of Ohlendorf et al. (2012) data
            % Output_pDusk(i,end-1)       = computed, normalized RFP_p readouts of ODE System
            % StDev_RFP_pDusk(i)          = standard deviation of Ohlendorf et al. (2012) data
        end
        S = sum(H(1:psetlength(1)));
    else
        for k = (psetlength(j-1)+1):psetlength(j);
            H = [H; ((Output_pDusk(k,3)-Output_pDusk(k,end-1))/Output_pDusk(k,end))^2];
            % time point is         = 17h (as always the end-point readout at 17h of the computation is store in Output
            % Output_pDusk(k,3)           = datapoints (RFP_p readout [FL/ OD600]) of Ohlendorf et al. (2012) data
            % Output_pDusk(k,end-1)       = computed, normalized RFP_p readouts of ODE System
            % StDev_RFP_pDusk(k)          = standard deviation of Ohlendorf et al. (2012) data
        end
        S = [S; sum(H((psetlength(j-1)+1):psetlength(j)))];
        % creating the scoring vector by summing up all individual scores per light intensity for each parameter set
    end    
end
% Column 1      = parameter_set_pDusk(j,:)
% Column 2      = N_pDusk
% Column 3      = RFP datapoints (pDusk)
% Column 4:10   = y1:y7
% Column 11     = normalized RFP_p (pDusk)
% Column 12     = StDev_RFP (pDusk)


%% Defining Score_pDusk
Score_pDusk = [1:j]';                                                           % scoring matrix with initial number of parameter sets
Score_pDusk = [Score_pDusk, S];                                                 % adding the S matrix from above


%%
%% pDawn
% Defining variables for pDawn
t   = 0:time_interval:17;   % Ohlendorf et al. (2013) induced the bacteria for 17h, i.e. 61200s. Thus the same time frame is used
x0  = [1 0 0 0 1 1 1 0 0];  % initial guesses for the variables of the ODE system. In the beginning state, the system is assumed to be in dark conditions, thus y_DD = 100% and j_a = 100% ("all" inactive FixJ (j_i) is phosphorylated)
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
N_pDawn = lightintensities_pDawn*10^-2;                 % �W*cm^-2 --- *10^-2 ---> W*m^-2
E_p = h*c/lambda_nm;                                    % distinct energy quanta of a photon
N_p = N_pDawn*lambda_nm*1E6/(1.988*1E-16);              % number of photons; lambda in [nm]
E_q = N_p/N_A;                                          % photonflux N = E_q in [�mol*m^-2*s^-1]
E_q_h = E_q*60*60;                                      % photonflux in [�mol*m^-2*h^-1]
N_pDawn = E_q_h;


%% Simulating Each parameter_set_pDawn(j) per Each Light Intensity N_pDawn(i)
X_Output = [];
P = [];
N_all = [];
RFP_all = [];
for j = 1:n
    for i = 1:length(N_pDawn)
        [T, X] = ode45(@(t,x) pDawn_function(t,x,parameter_set_pDawn(j,:),N_pDawn(i)), t, x0);
        X_Output = [X_Output; X(end,:)];
    end 
    J(1:length(N_pDawn)) = j;
    P = [P; J'];
    N_all = [N_all; N_pDawn']; 
    RFP_all = [RFP_all; RFP_datapoints_pDawn'];
end
Output_pDawn = [P, N_all, RFP_all, X_Output];     
% Column 1      = parameter_set_pDawn(j,:)
% Column 2      = N_pDawn
% Column 3      = RFP datapoints (pDawn)
% Column 4:12   = y1:y9
% Column 13     = normalized RFP_p (pDawn)
% Column 14     = StDev_RFP (pDawn)


%% Normalizing the Generated RFP_p Output (pDawn) 
% Here the generated RFP_p concentration is normlized according to the highest light intensity readout of RFP in Ohlendorf et al. (2012)
% level (FL/ OD600) given in Ohlendorf et al. (2012) figure 2a.
est_RFP_level = [];             % estimated RFP level
for h = 1:length(Output_pDawn(:,1))      
    est_RFP_level = [est_RFP_level; (0.894105/Output_pDawn(K(h),12))*Output_pDawn(h,12)]; 
    % 12                        = y(9)      = RFP_p
    % 0.894105                  = RFP datapoint from Ohlendorf et al. (2012) at N~100 (highest light intensity point) �W*cm^-2
    % Output_pDawn(K(h),12))    = taking the computed value of RFP_p = y(9) = Output(:,12) at the highest light intensity ~N100 of each parameter set
    % Output_pDawn(h,12)        = taking the computed value from each light intensity of each parameter set
end            
Output_pDawn = [Output_pDawn, est_RFP_level];
% Column 13     = normalized RFP_p (pDawn)


%% Loading the Standard Deviation Values from Ohlendorf et al. (2012)
F = csvread('Standard_Deviation_Ohlendorf_pDawn.csv');
StDev_Ohlendorf_pDawn = [F(:,2)'];                                                            % this is StDev in logscale
StDev_Sorted_pDawn = StDev_Ohlendorf_pDawn([1,3,5,7,9,11,13,15]);                             % subsetting upper bound StDev values
StDev_Sorted_pDawn = [StDev_Sorted_pDawn; StDev_Ohlendorf_pDawn([2,4,6,8,10,12,14,16])];      % subsetting lower bound StDev values
% StDev_RFP_pDawn = (StDev_Sorted_pDawn(1,:)-StDev_Sorted_pDawn(2,:))/2;                      % extracting StDev values (not considering log scale)
StDev_RFP_pDawn = (exp(StDev_Sorted_pDawn(1,:))-exp(StDev_Sorted_pDawn(2,:)))/2;              % extracting StDev values (considering log scale)
%!question: I think I need to transform the StDev from a log scale to a normal scale?!


%% Evaluating the Parameter Sets - The Scoring Function
% In order to understand whether a parameter set describes the given data
% by Ohlendorf et al. (2012) can be done by computing a scoring function.
% The idea is similar to the sum of squared errors (SSEs) approach and can
% be followed in more detail in Raue et al. (2009) equation IV.

% Adding the StDev_RFP to the matrix "Output"
Output_pDawn = [Output_pDawn, repmat(StDev_RFP_pDawn,1,j)'];
% Column 14     = StDev_RFP_pDawn

% The Scoring Function
% The agreement of experimental data with the observables predicted by the model is measured by an objective function, commonly the weighted sum of squared residuals - Raue et al. (2009)
psetlength = [1:j]*length(N_pDawn);
H = [];
S = [];
for j = 1:n
    if j == 1
        for i = 1:psetlength(1)
            H = [H; ((Output_pDawn(i,3)-Output_pDawn(i,end-1))/StDev_RFP_pDawn(i))^2];
            % time point is               = 17h (as always the end-point readout at 17h of the computation is store in Output
            % Output_pDawn(i,3)           = datapoints (RFP_p readout [FL/ OD600]) of Ohlendorf et al. (2012) data
            % Output_pDawn(i,end-1)       = computed, normalized RFP_p readouts of ODE System
            % StDev_RFP_pDawn(i)          = standard deviation of Ohlendorf et al. (2012) data
        end
        S = sum(H(1:psetlength(1)));
    else
        for k = (psetlength(j-1)+1):psetlength(j)
            H = [H; ((Output_pDawn(k,3)-Output_pDawn(k,end-1))/Output_pDawn(k,end))^2];
            % time point is               = 17h (as always the end-point readout at 17h of the computation is store in Output
            % Output_pDawn(k,3)           = datapoints (RFP_p readout [FL/ OD600]) of Ohlendorf et al. (2012) data
            % Output_pDawn(k,end-1)       = computed, normalized RFP_p readouts of ODE System
            % StDev_RFP_pDawn(k)          = standard deviation of Ohlendorf et al. (2012) data
        end
        S = [S; sum(H((psetlength(j-1)+1):psetlength(j)))];
    end    
end
% Column 1      = parameter_set_pDawn(j,:)
% Column 2      = N_pDawn
% Column 3      = RFP datapoints (pDawn)
% Column 4:12   = y1:y9
% Column 13     = normalized RFP_p (pDawn)
% Column 14     = StDev_RFP (pDawn)


%% Defining Score pDawn
Score_pDawn = [1:j]';                                                           % scoring matrix with initial number of parameter set
Score_pDawn = [Score_pDawn, S];                                                 % adding the S matrix from above



%%
%% Common Analysis

%% Weighted Sum Scoring Function
%% Finding the BestScore for both pDusk and pDawn
W = [];
Weighted_Sums = [];
for i = 1:n
    W = (Score_pDusk(i,2) + Score_pDawn(i,2))/2;                                % we do not apply a weight factor as there is no reason to assume that one system is more "important"
    Weighted_Sums = [Weighted_Sums W];
end
Weighted_Sums = [[1:j]' Weighted_Sums'];                                         % adding a count column to the Weighted_Sums
Score = Weighted_Sums;

BestScore = Weighted_Sums(Weighted_Sums(:,2) == min(Weighted_Sums(:,2)),:);     % finding the best score


%% Control Plot pDusk_pDawn
figure(2)
% pDusk
RFP_p_computed_pDusk = Output_pDusk(Output_pDusk(:,1) == BestScore(1), end-1); % taking the normalized concentrations of RFP_p depending on the ligth intensities of the parameter set with the best Score_pDusk
subplot(1,2,1);
% Plotting the Ohlendorf DataPoints for pDusk
errorbar(N_pDusk, RFP_datapoints_pDusk, StDev_RFP_pDusk,'.k')
set(gca,'YScale','log')
hold on
% Plotting the Computed RFP_p Levels for pDusk
scatter(N_pDusk, RFP_p_computed_pDusk)
plot(N_pDusk, RFP_p_computed_pDusk)
clear ylim;
ylim([10^-3 10]);
hold off

% pDawn
RFP_p_computed_pDawn = Output_pDawn(Output_pDawn(:,1) == BestScore(1), end-1); % taking the normalized concentrations of RFP_p depending on the ligth intensities of the parameter set with the best Score_pDawn
subplot(1,2,2);
% Plotting the Ohlendorf DataPoints for pDawn
errorbar(N_pDawn, RFP_datapoints_pDawn, StDev_RFP_pDawn,'.k')
set(gca,'YScale','log')
hold on
% Plotting the Computed RFP_p Levels for pDawn
scatter(N_pDawn, RFP_p_computed_pDawn)
plot(N_pDawn, RFP_p_computed_pDawn)
clear ylim;
ylim([10^-3 10]);
hold off


%% Write Outputs to csv files
% Parameter Set
spec = '%d%s%.3f%s%s%sParameter_Set.csv'; 
str = sprintf(spec, n, Underscore, time_interval, Underscore, System, Underscore);
dlmwrite(str, parameter_set);

% Score
spec = '%d%s%.3f%s%s%sScore.csv'; 
str = sprintf(spec, n, Underscore, time_interval, Underscore, System, Underscore);
dlmwrite(str, Score);

% Output_pDusk
spec = '%d%s%.3f%s%s%sOutput_pDusk.csv'; 
str = sprintf(spec, n, Underscore, time_interval, Underscore, System, Underscore);
dlmwrite(str, Output_pDusk);

% Output_pDawn
spec = '%d%s%.3f%s%s%sOutput_pDawn.csv'; 
str = sprintf(spec, n, Underscore, time_interval, Underscore, System, Underscore);
dlmwrite(str, Output_pDawn);

% Best_Output_pDusk
spec = '%d%s%.3f%s%s%sBestOutput_pDusk.csv'; 
str = sprintf(spec, n, Underscore, time_interval, Underscore, System, Underscore);
dlmwrite(str, RFP_p_computed_pDusk);

% Best_Output_pDawn
spec = '%d%s%.3f%s%s%sBestOutput_pDawn.csv'; 
str = sprintf(spec, n, Underscore, time_interval, Underscore, System, Underscore);
dlmwrite(str, RFP_p_computed_pDawn);

toc;