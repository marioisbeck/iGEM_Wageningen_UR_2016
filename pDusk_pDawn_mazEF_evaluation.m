%% pDusk_pDawn_mazEF_evaluation
%% Clear History
% clc;                % clear display
% clear all;          % clear Workspace
% close all;          % close all figures 
% clear workspace;    % clear workspace


%% 
%% Common Parameters of pDusk and pDawn
time_interval = 0.010;      % define time steps for the figures. 0.01h means taking a time point every 3.6s.
n_pDusk  = 23;              % number of parameters of the ODE system of pDusk
n_pDawn  = 28;              % number of parameters of the ODE system of pDawn
n_parameters = 28;          % n_parameters = n_pDawn (for better readability of code)
n = 1000;                   % number of the sample parameter sets to be read
System = 'pDusk_pDawn';     % for writing data files
Underscore = '_';           % for writing data files


%% Reading the Generated Data
Underscore = '_';

% Parameter Set pDusk
spec = '%d%s%.3f%s%s%sParameter_Set_pDusk.csv'; 
str = sprintf(spec, n, Underscore, time_interval, Underscore, System, Underscore);
parameter_set_pDusk = csvread(str);

% Parameter Set pDawn
spec = '%d%s%.3f%s%s%sParameter_Set_pDawn.csv'; 
str = sprintf(spec, n, Underscore, time_interval, Underscore, System, Underscore);
parameter_set_pDawn = csvread(str);

% Output_pDusk_mazEF
spec = '%d%s%.3f%s%s%sOutput_pDusk_mazEF.csv'; 
str = sprintf(spec, n, Underscore, time_interval, Underscore, System, Underscore);
Output_pDusk_mazEF = csvread(str);

% Output_pDawn_mazEF
spec = '%d%s%.3f%s%s%sOutput_pDawn_mazEF.csv'; 
str = sprintf(spec, n, Underscore, time_interval, Underscore, System, Underscore);
Output_pDawn_mazEF = csvread(str);

% pDusk Ratio
spec = '%d%s%.3f%s%s%spDusk_Ratio.csv'; 
str = sprintf(spec, n, Underscore, time_interval, Underscore, System, Underscore);
pDusk_Ratio = csvread(str);

% pDawn Ratio
spec = '%d%s%.3f%s%s%spDawn_Ratio.csv'; 
str = sprintf(spec, n, Underscore, time_interval, Underscore, System, Underscore);
pDawn_Ratio = csvread(str);


%% Constraining the A/T ratios 
% pDusk
constrained_pDusk_Ratio_bigger_one = pDusk_Ratio(pDusk_Ratio(:,3) > 1 & pDusk_Ratio(:,2) == 1,:);    % A/T > 1       at N(1)
constrained_pDusk_Ratio_smaller_one = pDusk_Ratio(pDusk_Ratio(:,3) < 1 & pDusk_Ratio(:,2) == 8,:);   % A/T < 1       at N(8)
subset = ismember(constrained_pDusk_Ratio_bigger_one(:,1), constrained_pDusk_Ratio_smaller_one(:,1));
subset = constrained_pDusk_Ratio_bigger_one(subset,1)';

constrained_pDusk_Ratio = [];
for i = 1:length(subset);
    constrained_pDusk_Ratio = [constrained_pDusk_Ratio; constrained_pDusk_Ratio_bigger_one(constrained_pDusk_Ratio_bigger_one(:,1) == subset(i),:)];
end
for i = 1:length(subset);
    constrained_pDusk_Ratio = [constrained_pDusk_Ratio; constrained_pDusk_Ratio_smaller_one(constrained_pDusk_Ratio_smaller_one(:,1) == subset(i),:)];
end

% pDawn
constrained_pDawn_Ratio_bigger_one = pDawn_Ratio(pDawn_Ratio(:,3) > 1 & pDawn_Ratio(:,2) == 1,:);    % A/T > 1       at N(1)
constrained_pDawn_Ratio_smaller_one = pDawn_Ratio(pDawn_Ratio(:,3) < 1 & pDawn_Ratio(:,2) == 8,:);   % A/T < 1       at N(8)

constrained_pDawn_Ratio = [];
for i = 1:length(subset);
    constrained_pDawn_Ratio = [constrained_pDawn_Ratio; constrained_pDawn_Ratio_bigger_one(constrained_pDawn_Ratio_bigger_one(:,1) == subset(i),:)];
end
for i = 1:length(subset);
    constrained_pDawn_Ratio = [constrained_pDawn_Ratio; constrained_pDawn_Ratio_smaller_one(constrained_pDawn_Ratio_smaller_one(:,1) == subset(i),:)];
end

% combining the constrains and storing the response
combined_constrains = [constrained_pDusk_Ratio, constrained_pDawn_Ratio(:,3)];


%% Storing the parameter sets per pDusk/ pDawn
p_pDusk = parameter_set_pDusk(subset(1:2),:);
p_pDawn = parameter_set_pDawn(subset(1:2),:);


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


%% Generating Response of mazE and mazF of pDusk per good parameter set
Xvalues = [];
Ni = [];
Time = [];
sub = [];
for k = 1:length(subset);
    for i = 1:length(N_pDusk);
        [T, X] = ode45(@(t,x) pDusk_function_const_mazF(t,x,p_pDusk(k,:),N_pDusk(i)), t, x0);
        Xvalues = [Xvalues; X];
    
        J = 1:length(X(:,1));
        J(J~=0) = i;
        Ni = [Ni; J'];
    
        Time = [Time; T];
        sub = [sub; repmat(subset(k), length(T), 1)];
    end 
end
X = [];
X_pDusk = [Ni, Time, sub, Xvalues(:,1:10)];


%% Plotting the response of mazE/ mazF of pDusk (subset(1))
figure(1);
for i = 1:length(N_pDusk);
    subplot(3,3,i);
	plot(T,X_pDusk(X_pDusk(:,1) == i & X_pDusk(:,3) == subset(1), [10 12]) );           % plotting mazE and mazF responses
    str = sprintf('N = %f',N_pDusk(i));
    str = {'Response of mazE and mazF in pDusk'; '+ constitutive mazF (ParameterSet 1)'; str};
    title(str);
    xlabel('time [h^1]');
    ylabel('concentration [µmol]');
end

% adding a legend to the graph
subplot(3,3,9);
white = plot(T(1),X_pDusk(1, [10 12]));
str = {'Legend'};
title(str);
legend ('e_p', 'f_p');                                                                          % plotting mazE and mazF responses


%% Plotting the response of mazE/ mazF of pDusk (subset(2))
figure(2);
for i = 1:length(N_pDusk);
    subplot(3,3,i);
    plot(T,X_pDusk(X_pDusk(:,1) == i & X_pDusk(:,3) == subset(2), [10 12]) );           % plotting mazE and mazF responses
    str = sprintf('N = %f',N_pDusk(i));
    str = {'Response of mazE and mazF in pDusk'; '+ constitutive mazF (ParameterSet 2)'; str};
    title(str);
    xlabel('time [h^1]');
    ylabel('concentration [µmol]');
end

% adding a legend to the graph
subplot(3,3,9);
white = plot(T(1),X_pDusk(1, [10 12]));
str = {'Legend'};
title(str);
legend ('e_p', 'f_p');                                                                          % plotting mazE and mazF responses


%%
%% pDawn
%% Loading the Data Points of Ohlendorf et al. (2012) for pDawn
% The data points of figure 2a in Ohlendorf et al. (2012) were retrieved
% through WebPlotDigitizer (http://arohatgi.info/WebPlotDigitizer) by (Ankit
% Rohatgi) to be used for the scoring function to evaluate the parameter
% sets generated by the latin hypercube sampling.
DataPoints_pDawn = csvread('DataPoints_pDawn.csv');     % reading the data points of figure 2a from Ohlendorf et al. (2012)
lightintensities_pDawn = [DataPoints_pDawn(:,1)'];      
RFP_datapoints_pDawn = [DataPoints_pDawn(:,2)'];

% Conversion of Light Intensities
N_pDawn = lightintensities_pDawn*10^-2;                 % µW*cm^-2 --- *10^-2 ---> W*m^-2
E_p = h*c/lambda_nm;                                    % distinct energy quanta of a photon
N_p = N_pDawn*lambda_nm*1E6/(1.988*1E-16);              % number of photons; lambda in [nm]
E_q = N_p/N_A;                                      % photonflux N = E_q in [µmol*m^-2*s^-1]
E_q_h = E_q*60*60;                                      % photonflux in [µmol*m^-2*h^-1]
N_pDawn = E_q_h;

% Defining variables for pDawn
t   = 0:time_interval:17;           % Ohlendorf et al. (2013) induced the bacteria for 17h, i.e. 61200s. Thus the same time frame is used
x0  = [1 0 0 0 1 1 1 0 0 0 0 0];    % initial guesses for the variables of the ODE system. In the beginning state, the system is assumed to be in dark conditions, thus y_DD = 100% and j_a = 100% ("all" inactive FixJ (j_i) is phosphorylated)
p = zeros(n_parameters);            % creating a zero-filled parameter vector


%% Generating Response of mazE and mazF of pDusk per good parameter set
Xvalues = [];
Ni = [];
Time = [];
sub = [];
for k = 1:length(subset);
    for i = 1:length(N_pDawn);
        [T, X] = ode45(@(t,x) pDawn_function_const_mazE(t,x,p_pDawn(k,:),N_pDawn(i)), t, x0);
        Xvalues = [Xvalues; X];
    
        J = 1:length(X(:,1));
        J(J~=0) = i;
        Ni = [Ni; J'];
    
        Time = [Time; T];
        sub = [sub; repmat(subset(k), length(T), 1)];
    end 
end
X = [];
X_pDawn = [Ni, Time, sub, Xvalues(:,1:12)];


%% Plotting the response of mazE/ mazF of pDawn (subset(1))
figure(3);
for i = 1:length(N_pDawn);
    subplot(3,3,i);
    plot(T,X_pDawn(X_pDawn(:,1) == i & X_pDawn(:,3) == subset(1), [14 12]) );           % plotting mazE and mazF responses
    str = sprintf('N = %f',N_pDusk(i));
    str = {'Response of mazE and mazF in pDawn'; '+ constitutive mazE (ParameterSet 1)'; str};
    title(str);
    xlabel('time [h^1]');
    ylabel('concentration [µmol]');
    clear ylim;
    ylim([0 (4*10^(-3))]);
end

% adding a legend to the graph
subplot(3,3,9);
white = plot(T(1),X_pDawn(1, [14 12]));
str = {'Legend'};
title(str);
legend ('e_p', 'f_p');                                                                          % plotting mazE and mazF responses
clear ylim;
ylim([0 (4*10^(-3))]);


%% Plotting the response of mazE/ mazF of pDawn (subset(2))
figure(4);
for i = 1:length(N_pDawn);
    subplot(3,3,i);
    plot(T,X_pDawn(X_pDawn(:,1) == i & X_pDawn(:,3) == subset(2), [14 12]) );           % plotting mazE and mazF responses
    str = sprintf('N = %f',N_pDusk(i));
    str = {'Response of mazE and mazF in pDawn'; '+ constitutive mazE (ParameterSet 2)'; str};
    title(str, 'FontSize', 14);
    xlabel('Time [h^1]', 'FontSize', 12);
    ylabel('Concentration [µmol]', 'FontSize', 12);
    clear ylim;
    ylim([0 0.03]);
    set(findall(gca, 'Type', 'Line'),'LineWidth',1);
    res = 300;
end

% adding a legend to the graph
subplot(3,3,9);
white = plot(T(1),X_pDawn(1, [13 11]));
str = {'Legend'};
title(str);
legend ('e_p', 'f_p');                                                                                          % plotting mazE and mazF responses
clear ylim;
ylim([0 0.03]);
set(gca, 'XTick', [], 'YTick', [])
set(findall(gca, 'Type', 'Line'),'LineWidth',1);
res = 300;

print('-painters', '-dpdf', '-f4');
%% T--Wageningen_UR--pDusk_pDawn_AT_Response_1
% Plotting the response of mazE/ mazF of pDusk (subset(1))
for i = 1:8;
figure(i);
subplot(1,2,1);
plot(T,X_pDusk(X_pDusk(:,1) == i & X_pDusk(:,3) == subset(1), [10 12]) );           % plotting mazE and mazF responses
str = sprintf('N = %f',N_pDusk(i));
str = {'pDusk + constitutive mazF'; str};
title(str, 'FontSize', 14);
xlabel('Time [h^1]', 'FontSize', 12);
ylabel('Concentration [µmol]', 'FontSize', 12);
clear ylim;
ylim([0 0.03]);
hold on

plot(T,X_pDusk(X_pDusk(:,1) == i & X_pDusk(:,3) == subset(2), [10 12]) );           % plotting mazE and mazF responses
title(str, 'FontSize', 14);
xlabel('Time [h^1]', 'FontSize', 12);
ylabel('Concentration [µmol]', 'FontSize', 12);
clear ylim;
ylim([0 0.03]);
legend ('e_p - set_1', 'f_p - set_1', 'e_p - set_2', 'f_p - set_2');
hold off

subplot(1,2,2);
plot(T,X_pDawn(X_pDawn(:,1) == i & X_pDawn(:,3) == subset(1), [14 12]) );           % plotting mazE and mazF responses
str = sprintf('N = %f.2',N_pDawn(i));
str = {'pDawn + constitutive mazE'; str};
title(str, 'FontSize', 14);
xlabel('Time [h^1]', 'FontSize', 12);
ylabel('Concentration [µmol]', 'FontSize', 12);
clear ylim;
ylim([0 0.03]);
hold on

plot(T,X_pDawn(X_pDawn(:,1) == i & X_pDawn(:,3) == subset(2), [14 12]) );           % plotting mazE and mazF responses
title(str, 'FontSize', 14);
xlabel('Time [h^1]', 'FontSize', 12);
ylabel('Concentration [µmol]', 'FontSize', 12);
clear ylim;
ylim([0 0.03]);
set(findall(gca, 'Type', 'Line'),'LineWidth',1);
res = 300;
hold off
end

print('-painters', '-dpdf', '-f1');
print('-painters', '-dpdf', '-f2');
print('-painters', '-dpdf', '-f3');
print('-painters', '-dpdf', '-f4');
print('-painters', '-dpdf', '-f5');
print('-painters', '-dpdf', '-f6');
print('-painters', '-dpdf', '-f7');
print('-painters', '-dpdf', '-f8');
%% T--Wageningen_UR--pDusk_pDawn_AT_Ratio
% Plotting the A/T Ratios of pDusk and pDawn
% pDusk E = 10/ F = 12
% pDawn E = 14/ F = 12
figure(5);
% pDusk
subplot(1,2,1);
for k = 1:length(subset);
    plot(N_pDusk, pDusk_Ratio(pDusk_Ratio(:,1) == subset(k), 3));
	set(gca,'YScale','log');
    clear ylim;
    ylim([10^-2 10^3]);
    set(findall(gca, 'Type', 'Line'),'LineWidth',1);
    res = 300;
    hold on
end
str = ('per Increasing Light Intensity N');
str = {'Decrease of pDusk A/T Ratios'; str};
title(str, 'FontSize', 14);
xlabel('Light Intensities N [µmol*m^-2*h^-1]', 'FontSize', 12);
ylabel('mazE/ mazF = Antitoxin/ Toxin (A/T)', 'FontSize', 12);
set(findall(gca, 'Type', 'Line'),'LineWidth',1);
res = 300;
hold off

% pDawn
subplot(1,2,2);
for k = 1:length(subset);
    plot(N_pDawn, pDawn_Ratio(pDawn_Ratio(:,1) == subset(k), 3));
    set(gca,'YScale','log');
    hold on
end
str = ('per Increasing Light Intensity N');
str = {'Decrease of pDusk A/T Ratios'; str};
title(str, 'FontSize', 14);
xlabel('Light Intensities N [µmol*m^-2*h^-1]', 'FontSize', 12);
ylabel('mazE/ mazF = Antitoxin/ Toxin (A/T)', 'FontSize', 12);
set(findall(gca, 'Type', 'Line'),'LineWidth',1);
res = 300;
hold off

print('-painters', '-dpng', '-f5');