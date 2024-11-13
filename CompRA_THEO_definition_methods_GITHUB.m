clear all
close all
clc
%startTime = tic;
startTimeCPU = cputime;

% Parameterize the life cycles of 5 distinct species.

Data = {[11,0.0365,0.0965,0.9],[5, 0.2, 0.25, 0.572],[1,0.8,0.385,0.4],[0.3, 0.93, 0.505, 0.3], [0.232, 0.95, 0.80, 0.07]}; %  (F, SA, SJ, Y)
nspecies = length(Data); % number of species

% Quantitative genetics information
Vp = 1; % Phenotypic variance
h2 = 0.2; % Heritability value of the phenotypic trait
Va = h2*Vp; % Additive genetic variance (h2 = Va/Vp) VA : 0.05, 0.2, 0.8
Ve = Vp - Va; % Environmental variance (Vp = Va + Ve)

% Set selective pressure (slope of the relationship between vital rate and phenotype)
beta = 0.15; %B : 0.01, 0.15, 0.5 

% Get demographic descriptive info
DEMOout=zeros(nspecies, 3); % outputs lambda, R0, Tgen
RAout=zeros(nspecies, 24); % outputs RA per unit of time, generation ,selection theta & R0)
S=zeros(nspecies,5); % sensitivity d lambda/ d theta
ndemo = 24;
for ispecies=1:nspecies
    % Extract vital rates of the species
    F  = Data{ispecies}(1);% This is the total number of offspring produced by a female per year that survive the 1st year
    SA = Data{ispecies}(2);
    SJ = Data{ispecies}(3);
    Y  = Data{ispecies}(4);
    f  = F/SJ; % This is the total number of offspring produced by a female per year

    % Projection matrix
    A = [SJ*(1-Y), f*SJ; SJ*Y, SA];

    % pop growth rate / fitness
    lambda = (SJ*(1-Y)+SA)/2+sqrt((SJ*(1-Y)-SA).^2+4*f*(SJ^2)*Y)/2;

    % eigenvectors to calculate sensitivty
    v = [Y*SJ lambda-SJ*(1-Y)]';
    w = [f*SJ lambda-SJ*(1-Y)]';

    % Get descriptive demographic information

    % life time reproductive success
    R0     = Y.*f*(SJ^2)./((1-SA)*(1-SJ*(1-Y)));
    % generation time
    Tgen = lambda*sum(v.*w) / (v(1)*w(2)*f*SJ);
    % D lambda / D A
    Sensitivity= v*w'./sum(v.*w);
    % Sensitivity lower level derivative d A/ d theta
    theta=[f SJ SA Y];
    PDlv=zeros(2,2,4);
    PDlv(1,2,1)=theta(2); % d A/ d theta wrt f
    PDlv(1,1,2)=(1-theta(4)); % d A/ d theta wrt SJ
    PDlv(1,2,2)=theta(1);% d A/ d theta wrt SJ
    PDlv(2,1,2)=theta(4);% d A/ d theta wrt SJ
    PDlv(2,2,3)=1; % d A/ d theta wrt SA
    PDlv(1,1,4)=-theta(2); % d A/ d theta wrt Y
    PDlv(2,1,4)=theta(2);% d A/ d theta wrt Y

    for j=1:length(theta)
        S(ispecies,j)=sum(sum(PDlv(:,:,j).*Sensitivity));
    end

    % calculate sensitivity of lambda to F
    temp=zeros(2,2); temp(1,2)=1;
    S(ispecies, 5)=sum(sum(temp.*Sensitivity));


    % STORE OUTPUTS in RAout
    % demographic paramters
    DEMOout(ispecies, 1) = lambda;
    DEMOout(ispecies, 2) = R0;
    DEMOout(ispecies, 3) = Tgen;
    % adaptation rates ********************************
    % --------------------------------------------------
    % fecundity
    % --------------------------------------------------
    % per time
    RAout(ispecies, 1) = h2*Vp*beta/ Tgen; %vi f time
    RAout(ispecies, 2) = h2*Vp/ Tgen; %vi f time unit selection theta
    RAout(ispecies, 3) = h2*Vp/ Tgen; %vi f time unit selection R0
    % per generation
    RAout(ispecies, 4) = h2*Vp*beta; %vi f gen
    RAout(ispecies, 5) = h2*Vp; %vi f gen unit selection theta
    RAout(ispecies, 6) = h2*Vp; %vi f gen unit selection fitness

    % --------------------------------------------------
    % juvenile survival
    % --------------------------------------------------
    % per time
    RAout(ispecies, 7) = h2*Vp*beta*(1/ Tgen)*(1+(lambda/(lambda-SJ*(1-Y))))*(1-SJ); %vi SJ time
    RAout(ispecies, 8) = h2*Vp*(1/Tgen)*(1+(lambda/(lambda-SJ*(1-Y)))); %vi SJ time unit selection theta
    RAout(ispecies, 9) = h2*Vp*(1/Tgen)*((1+(lambda/(lambda-SJ*(1-Y))))/(1+(1/(1-SJ*(1-Y))))); %vi SJ time unit selection R0
    % per generation
    RAout(ispecies, 10) = h2*Vp*beta*(1+(lambda/(lambda-SJ*(1-Y))))*(1-SJ); %vi SJ gen
    RAout(ispecies, 11) = h2*Vp*(1+(lambda/(lambda-SJ*(1-Y)))); %vi SJ gen unit selection theta
    RAout(ispecies, 12) = h2*Vp*((1+(lambda/(lambda-SJ*(1-Y))))/(1+(1/(1-SJ*(1-Y))))); %vi SJ time unit selection R0


    % --------------------------------------------------
    % maturation rate
    % --------------------------------------------------
    % per time
    RAout(ispecies, 13) =h2*Vp*beta*(1/Tgen)*(((lambda-SJ)*(1-Y)))/(lambda-(SJ*(1-Y))); %vi Y time
    RAout(ispecies, 14) =h2*Vp*(1/Tgen)*(lambda-SJ)/(lambda-SJ*(1-Y));%*(1-(Y*lambda/(lambda-SJ*(1-Y)))); %vi Y time unit selection theta
    RAout(ispecies, 15) =h2*Vp*(1/Tgen)*((lambda-SJ)*(1-SJ*(1-Y))/((1-SJ)*(lambda-SJ*(1-Y))));%vi Y time unit selection R0
    % per generation
    RAout(ispecies, 16) =h2*Vp*beta*(((lambda-SJ)*(1-Y)))/(lambda-(SJ*(1-Y))); %vi Y gen
    RAout(ispecies, 17) =h2*Vp*(lambda-SJ)/(lambda-SJ*(1-Y));%(1- (Y*lambda/(lambda-SJ*(1-Y))));%vi Y gen unit selection theta
    RAout(ispecies, 18) =h2*Vp*((lambda-SJ)*(1-SJ*(1-Y))/((1-SJ)*(lambda-SJ*(1-Y))));%vi Y gen unit selection R0


    % --------------------------------------------------
    % adult survival
    % --------------------------------------------------
    % per time
    RAout(ispecies, 19) = h2*Vp*beta*(1/Tgen)* ((SA*(1-SA))/((lambda-SA)));  % vi Sa time
    RAout(ispecies, 20) = h2*Vp*(1/Tgen)*(SA/(lambda-SA));  % vi Sa time unit selection theta
    RAout(ispecies, 21) = h2*Vp*(1/Tgen)*((1-SA)/(lambda-SA));  % vi Sa time unit selection R0
    % per generation
    RAout(ispecies, 22) = h2*Vp*beta*((SA*(1-SA))/((lambda-SA))); %vi Sa gen
    RAout(ispecies, 23) = h2*Vp*(SA/(lambda-SA)); %vi Sa gen unit selection theta
    RAout(ispecies, 24) = h2*Vp*((1-SA)/(lambda-SA)); %vi Sa gen unit selection R0

end

endTimeCPU = cputime - startTimeCPU;  % Calculate CPU elapsed time
fprintf('CPU Time: %f seconds\n', endTimeCPU);

%%
100*(RAout(:, [1 7 19 13]))

%% FIGURE WITH 6 DEFINITIONS
species_labels = {'species 1', 'species 2', 'species 3', 'species 4', 'species 5'};

colors_rgb = [
    255, 192, 203;   % Pink
    0, 100, 0;       % Dark Green
    0, 0, 255;       % Blue
    128, 0, 128;     % Purple
    255, 165, 0      % Orange
    ] / 255; % Scale RGB values to [0, 1]

% Plotting

figure1 = figure('WindowState','fullscreen',...
    'Colormap',[1 0.75 0.79;0 0.39 0;0 0 1;0.50 0 0.50;1 0.64 0]);

subplot(3,2,1)
% Create multiple lines using matrix input to bar
bar1 = bar(RAout(:, 1:6:ndemo));
set(bar1(1),'FaceColor',[1 0 1]);
set(bar1(2),...
    'FaceColor',[0.49 0.18 0.55]);
set(bar1(4),...
    'FaceColor',[0.30 0.74 0.93]);

legend('f','SJ','Y','SA')
%xlabel('Vital Rates');
ylabel('RA_{T}');
title('Rate of adaptation per time');
xticks(1:5); % Adjust ispecies-axis ticks to match the number of vital rates
xticklabels(species_labels); % Set ispecies-axis tick labels
set(gca, 'FontSize', 24);

subplot(3,2,2)
% Create multiple lines using matrix input to bar
bar1 = bar(RAout(:, 4:6:ndemo));
set(bar1(1),'FaceColor',[1 0 1]);
set(bar1(2),...
    'FaceColor',[0.49 0.18 0.55]);
set(bar1(4),...
    'FaceColor',[0.30 0.74 0.93]);

%legend('f','SJ','Y','SA')
%xlabel('Vital Rates');
ylabel('RA_{G}');
title('Rate of adaptation per generation');
xticks(1:5); % Adjust ispecies-axis ticks to match the number of vital rates
xticklabels(species_labels); % Set ispecies-axis tick labels
set(gca, 'FontSize', 24);

subplot(3,2,3)
% Create multiple lines using matrix input to bar
bar1 = bar(RAout(:, 2:6:ndemo));
set(bar1(1),'FaceColor',[1 0 1]);
set(bar1(2),...
    'FaceColor',[0.49 0.18 0.55]);
set(bar1(4),...
    'FaceColor',[0.30 0.74 0.93]);

%legend('f','SJ','Y','SA')
%xlabel('Vital Rates');
ylabel('RA_{TS_\theta}');
title('Rate of adaptation per time per unit of selection THETA');
xticks(1:5); % Adjust ispecies-axis ticks to match the number of vital rates
xticklabels(species_labels); % Set ispecies-axis tick labels
set(gca, 'FontSize', 24);

subplot(3,2,4)
% Create multiple lines using matrix input to bar
bar1 = bar(RAout(:, 5:6:ndemo));
set(bar1(1),'FaceColor',[1 0 1]);
set(bar1(2),...
    'FaceColor',[0.49 0.18 0.55]);
set(bar1(4),...
    'FaceColor',[0.30 0.74 0.93]);

%legend('f','SJ','Y','SA')
%xlabel('Vital Rates');
ylabel('RA_{GS_\theta}');
title('Rate of adaptation per generation per unit of selection THETA');
xticks(1:5); % Adjust ispecies-axis ticks to match the number of vital rates
xticklabels(species_labels); % Set ispecies-axis tick labels
set(gca, 'FontSize', 24);

subplot(3,2,5)
% Create multiple lines using matrix input to bar
bar1 = bar(RAout(:, 3:6:ndemo));
set(bar1(1),'FaceColor',[1 0 1]);
set(bar1(2),...
    'FaceColor',[0.49 0.18 0.55]);
set(bar1(4),...
    'FaceColor',[0.30 0.74 0.93]);

%legend('f','SJ','Y','SA')
%xlabel('Vital Rates');
ylabel('RA_{TSR_0}');
title('Rate of adaptation per time per unit of selection R0');
xticks(1:5); % Adjust ispecies-axis ticks to match the number of vital rates
xticklabels(species_labels); % Set ispecies-axis tick labels
set(gca, 'FontSize', 24);

subplot(3,2,6)
% Create multiple lines using matrix input to bar
bar1 = bar(RAout(:, 6:6:ndemo));
set(bar1(1),'FaceColor',[1 0 1]);
set(bar1(2),...
    'FaceColor',[0.49 0.18 0.55]);
set(bar1(4),...
    'FaceColor',[0.30 0.74 0.93]);

%legend('f','SJ','Y','SA')
%xlabel('Vital Rates');
ylabel('RA_{GSR_0}');
title('Rate of adaptation per generation per unit of selection R0');
xticks(1:5); % Adjust ispecies-axis ticks to match the number of vital rates
xticklabels(species_labels); % Set ispecies-axis tick labels
set(gca, 'FontSize', 24);
%%  comparaison methods
species_labels = {'species 1', 'species 2', 'species 3', 'species 4', 'species 5'};
vital_rate_labels = {'f THEO',  'f IBM','f MPM',};

figure(2)
subplot(2,2,1)

demo_dataf=[1.47	1.55	1.45
1.27	1.28	1.27
0.48	0.51	0.49
0.19	0.27	0.20
0.13	0.16	0.14];

% Create multiple lines using matrix input to bar
bar1 = bar(demo_dataf');
xlabel('Vital Rates');
ylabel('fecundity');
title('Fecundity');
xticks(1:3); % Adjust ispecies-axis ticks to match the number of vital rates
xticklabels(vital_rate_labels); % Set ispecies-axis tick labels
set(gca, 'FontSize', 24);
ylim([0 1.6])

subplot(2,2,2)
demo_dataSJ=[2.662	2.422	2.4217
    2.013	1.874	1.8736
    0.674	0.681	0.6808
    0.240	0.246	0.2458
    0.124	0.128	0.1282];

vital_rate_labels = {'SJ THEO','SJ IBM', 'SJ MPM'};
bar(demo_dataSJ');
xlabel('Vital Rates');
ylabel('juvenile survival');
title('Juvenile survival');
xticks(1:3); % Adjust ispecies-axis ticks to match the number of vital rates
xticklabels(vital_rate_labels); % Set ispecies-axis tick labels
set(gca, 'FontSize', 24);

subplot(2,2,3)
demo_dataY=[0.134	0.131	0.1308
    0.455	0.436	0.436
    0.229	0.221	0.2213
    0.102	0.1	0.1002
    0.092	0.088	0.0885];

vital_rate_labels = {'Y THEO', 'Y IBM','Y MPM'};
bar(demo_dataY');
xlabel('Vital Rates');
ylabel('maturation');
title('Maturation rate');
xticks(1:3); % Adjust ispecies-axis ticks to match the number of vital rates
xticklabels(vital_rate_labels); % Set ispecies-axis tick labels
set(gca, 'FontSize', 24);

subplot(2,2,4)
demo_dataSA=[0.053	0.053	0.0529
    0.253	0.255	0.2545
    0.381	0.362	0.3616
    0.176	0.168	0.1686
    0.119	0.113	0.1131];

vital_rate_labels = {'SA THEO', 'SA IBM', 'SA MPM'};
bar(demo_dataSA');
xlabel('Vital Rates');
ylabel('adult survival');
title('Adult survival');
xticks(1:3); % Adjust ispecies-axis ticks to match the number of vital rates
xticklabels(vital_rate_labels); % Set ispecies-axis tick labels
set(gca, 'FontSize', 24);
legend(species_labels, 'FontSize', 24)


%%
% Plotting like R
fo=24;

% Define colors for vital rates using ggplot2 palette
colors = [
    248/255, 118/255, 109/255; % 'f' - fertility (reddish)
    124/255, 174/255, 0/255;   % 'SJ' - juvenile survival (greenish)
    0/255,   191/255, 196/255; % 'Y' - maturation (bluish)
    199/255, 124/255, 255/255; % 'SA' - adult survival (purple)
];

% Open figure in fullscreen mode
figure1 = figure('WindowState','fullscreen');

% Assuming you have species labels defined
% Replace with your actual species labels
%species_labels = {'Species 1', 'Species 2', 'Species 3', 'Species 4', 'Species 5'};
species_labels = {'1', '2', '3', '4', '5'};

% Adjust ndemo if not already defined
% ndemo = total number of columns in RAout
% For example, if RAout has 30 columns:
% ndemo = 30;

% Create legend entries only once
legend_entries = {'f', 'SJ', '\gamma', 'SA'};

% Subplot 1
subplot(3,2,1)
% Create bar plot
bar1 = bar(RAout(:, 1:6:ndemo), 'grouped');
% Set colors for each vital rate
for idx = 1:length(bar1)
    set(bar1(idx), 'FaceColor', colors(idx,:));
end
% Set labels and title
ylabel('RA_{T}', 'FontSize', fo);
title('Rate of adaptation per time', 'FontSize', fo);
xticks(1:5);
xticklabels(species_labels);
set(gca, 'FontSize', fo);
% Legend
legend_handle = legend(legend_entries, 'Location', 'eastoutside', 'FontSize', fo);
grid on; % Add grid
% Subplot 2
subplot(3,2,2)
bar1 = bar(RAout(:, 4:6:ndemo), 'grouped');
for idx = 1:length(bar1)
    set(bar1(idx), 'FaceColor', colors(idx,:));
end
ylabel('RA_{G}', 'FontSize', fo);
title('Rate of adaptation per generation', 'FontSize', fo);
xticks(1:5);
xticklabels(species_labels);
set(gca, 'FontSize', fo);
grid on; % Add grid
% Subplot 3
subplot(3,2,3)
bar1 = bar(RAout(:, 2:6:ndemo), 'grouped');
for idx = 1:length(bar1)
    set(bar1(idx), 'FaceColor', colors(idx,:));
end
ylabel('RA_{TS_\theta}', 'FontSize', fo);
title('Rate of adaptation per time per unit of selection \theta', 'FontSize', fo);
xticks(1:5);
xticklabels(species_labels);
set(gca, 'FontSize', fo);
grid on; % Add grid
% Subplot 4
subplot(3,2,4)
bar1 = bar(RAout(:, 5:6:ndemo), 'grouped');
for idx = 1:length(bar1)
    set(bar1(idx), 'FaceColor', colors(idx,:));
end
ylabel('RA_{GS_\theta}', 'FontSize', fo);
title('Rate of adaptation per generation per unit of selection \theta', 'FontSize', fo);
xticks(1:5);
xticklabels(species_labels);
set(gca, 'FontSize', fo);
grid on; % Add grid
% Subplot 5
subplot(3,2,5)
bar1 = bar(RAout(:, 3:6:ndemo), 'grouped');
for idx = 1:length(bar1)
    set(bar1(idx), 'FaceColor', colors(idx,:));
end
ylabel('RA_{TSR_0}', 'FontSize', fo);
title('Rate of adaptation per time per unit of selection R_0', 'FontSize', fo);
xticks(1:5);
xticklabels(species_labels);
set(gca, 'FontSize', fo);
xlabel('Species')
grid on; % Add grid
% Subplot 6
subplot(3,2,6)
bar1 = bar(RAout(:, 6:6:ndemo), 'grouped');
for idx = 1:length(bar1)
    set(bar1(idx), 'FaceColor', colors(idx,:));
end
ylabel('RA_{GSR_0}', 'FontSize', fo);
title('Rate of adaptation per generation per unit of selection R_0', 'FontSize', fo);
xticks(1:5);
xticklabels(species_labels);
set(gca, 'FontSize', fo);
grid on; % Add grid
% Add a super title for the figure
sgtitle('Comparison of adaptation rates based on six definitions', 'FontSize', fo);
xlabel('Species')
% Adjust the position of the legend
legend_handle.Position = [0.92, 0.4, 0.05, 0.2]; % Adjust as needed

% Adjust subplots to make room for the legend
for i = 1:6
    subplot(3,2,i);
    pos = get(gca, 'Position');
    pos(3) = pos(3) - 0.05; % Reduce width to make room for legend
    set(gca, 'Position', pos);
end