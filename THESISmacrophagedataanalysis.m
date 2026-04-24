%this is where i will attempt to analyze the macrophage population data
%% UNTREATED CONTROL TRIALS
clear all
load('kyanna30trialsdata.mat'); %this data is for the 30 trials data, representing the untreated controls

untreatedcontrol_NP = nSTORAGE(:,:);
untreatedcontrol_M1 = m1STORAGE(:,:);
untreatedcontrol_M2 = m2STORAGE(:,:);
%% 
close all

%Plotting M1 population vs time
%time = linspace(0, 1001, 1001); %this is if I want time to be in terms of timesteps (data saved every 60 seconds of the 60,000 minute simulation, makes for 1000 time steps
time = (0:1000) * 60 / 1440; %this converts time to simulated 41.67 days for 60,000 minutes
figure;

tl = tiledlayout(1,2); % 1 row, 2 columns
nexttile(tl, 1)  %left tile

title("Untreated Control M1 Population vs Time")
xlabel("Simulation Time (days)")
ylabel("Cell Count")
hold on; 
for i = 4:33
    plot(time, untreatedcontrol_M1(i, :), "Color",[0.7 0.6 0.9]);
end
xlim([0 41.67])



%data analysis
untreatedcontrol_M1_mean_curve = mean(untreatedcontrol_M1, 1); %1 means it is averaging across rows 
plot(time, untreatedcontrol_M1_mean_curve, 'LineWidth', 2, 'Color', [0 0 0]);

%bootstrapping a 95% confidence interval of the mean
untreatedcontrol_M1_data_to_boot = untreatedcontrol_M1(4:33, :); %same size we want to bootstrap
untreatedcontrol_M1_ci = bootci(1000,@mean,untreatedcontrol_M1_data_to_boot); 
plot(time, untreatedcontrol_M1_ci(1,:), 'LineWidth', 0.75, 'Color', [0 0 0], 'LineStyle','--');
plot(time, untreatedcontrol_M1_ci(2,:), 'LineWidth', 0.75, 'Color', [0 0 0], 'LineStyle','--');

%95% CI upper and lower bounds
ci_lower = untreatedcontrol_M1_ci(1,:);
ci_upper = untreatedcontrol_M1_ci(2,:);

%shading confidence interval
fill([time fliplr(time)], ...
     [ci_upper fliplr(ci_lower)], ...
     'y', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

%Plotting M2
time = linspace(0, 1001, 1001);
nexttile(tl,2)

title("Untreated Control M2 Population vs Time")
xlabel("Simulation Time (days)")
ylabel("Cell Count")
hold on; 
for i = 4:33
    plot(time, untreatedcontrol_M2(i, :), "Color",[0.7 0.6 0.9]);
end
xlim([0 41.67])

%bootstrapping w/ 95% CI to see what the differnce is
untreatedcontrol_M2_mean_curve = mean(untreatedcontrol_M2, 1); %1 means it is averaging across rows 
plot(time, untreatedcontrol_M2_mean_curve, 'LineWidth', 2, 'Color', [0 0 0]);

%bootstrapping a 95% confidence interval of the mean
untreatedcontrol_M2_data_to_boot = untreatedcontrol_M2(4:33, :); %same size we want to bootstrap
untreatedcontrol_M2_ci = bootci(1000,@mean,untreatedcontrol_M2_data_to_boot); 
plot(time, untreatedcontrol_M2_ci(1,:), 'LineWidth', 0.75, 'Color', [0 0 0], 'LineStyle','--');
plot(time, untreatedcontrol_M2_ci(2,:), 'LineWidth', 0.75, 'Color', [0 0 0], 'LineStyle','--');

% 95% CI upper and lower bounds (for shading purposes)
ci_lower = untreatedcontrol_M2_ci(1,:);
ci_upper = untreatedcontrol_M2_ci(2,:);

%shading confidence interval
fill([time fliplr(time)], ...
     [ci_upper fliplr(ci_lower)], ...
     'y', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

%Comparing M1 and M2 means on same graph
figure;
plot(time, untreatedcontrol_M1_mean_curve, 'LineWidth', 2, 'Color', 'r');
hold on
plot(time, untreatedcontrol_M2_mean_curve, 'LineWidth', 2, 'Color', 'b');
title("Average M1 and M2 Populations Over Time - Untreated Control")
xlabel("Simulation Time (days)")
ylabel("Cell Count")
legend('M1','M2')
%% trying to plot m1 and m2 for each trial independently to see if we can see any oscilations
%and also counting number of change overs
intersection_counts = zeros(1, 30); %  zero vectorstore counts for each trial
all_intervals = []; %zero vector

%USING TILEDLAYOUT
figure;
tl = tiledlayout(6,5);
title(tl,'Untreated Control M1 and M2 Populations vs Time');
xlabel(tl,'Simulation Time (days)');
ylabel(tl,'Cell Count');


for i = 4:33
   
    nexttile(tl)
    plot(time, untreatedcontrol_M1(i, :), "Color",'r');
    hold on;
    plot(time, untreatedcontrol_M2(i, :), "Color",'b');

    title(['Trial', num2str(i-3)]);
    
    xlim([0 41.67])


    %Intersection stuff
    diff_signal = untreatedcontrol_M1(i,:) - untreatedcontrol_M2(i,:);

    % Find where the sign changes
    sign_changes = diff_signal(1:end-1) .* diff_signal(2:end) < 0;

    intersection_counts(i-3) = sum(sign_changes);
   
    %Trying to determine the time intervals between intersections
    % Convert indices to times
    t_cross = time(sign_changes);

    % Compute time intervals between crossings
    dt = diff(t_cross);

    % Store them
    all_intervals = [all_intervals, dt];

     %existing intersection count
    count = intersection_counts(i-3);

    % Add text annotation inside the plot
text(0.98, 0.98, ['Dominance shifts: ', num2str(count)], ...
    'Units', 'normalized', ...         % position relative to axes (0-1)
    'HorizontalAlignment', 'right', ...% align text to the right
    'VerticalAlignment', 'top', ...    % align text to the top
    'FontWeight', 'bold', ...
    'FontSize', 10, ...
    'Margin', 2);                      % small padding inside box

end
legend('M1', 'M2','Location','eastoutside');

% TRYING TO PLOT M2:M1 RATIO OVER TIME FOR ALL TRIALS
ratios_overtime = []; %initializing zero vector for thing i want to calculate/plot

%USING TILEDLAYOUT
figure;
tl = tiledlayout(6,5);
title(tl,'M2:M1 Ratio in Untreated Control vs Time');
xlabel(tl,'Simulation Time (days)');
ylabel(tl,'M2:M1 Ratio');
dt = 60;
time_days = (0:1000) * dt / 1440;

for i = 4:33
  
    nexttile(tl)
    %calculate ratio
    ratios_overtime = untreatedcontrol_M2(i, :)./(untreatedcontrol_M1(i, :) + untreatedcontrol_M2(i, :)) ;
    plot(time_days, ratios_overtime, "Color",[0.231 0.208 0.380]);
  

    title(['Trial ', num2str(i-3)]);
    
    xlim([0 42])
    ylim([0.3 0.7])

    all_ratios(i-3,:) = ratios_overtime;
end

%making plot of mean curve across all trials
figure;
untreated_m2ratio_mean_curve = mean(all_ratios,1);
plot(time_days, untreated_m2ratio_mean_curve,'k','LineWidth',2,"Color",[0.231 0.208 0.380]);
xlabel("Time (Simulation days)")
ylabel("M2:M1 Ratio")
title("Mean M2:M1 Ratio vs Simulated Time")

%making a histogram of the number of intersections
figure;
histogram(intersection_counts);

title('Distribution of Intersection Counts Across Untreated Control Trials');
xlabel('Number of Intersections');
ylabel('Frequency');

%making histogram of timing betwe osciallations
figure;
histogram(all_intervals);

title('Distribution of Time Between Intersections - Untreated Control');
xlabel('Time Between Crossings');
ylabel('Frequency');

%measure how regular the oscillations are using the coefficient of variation (CV) of the intervals:

CV = std(all_intervals) / mean(all_intervals);

% Relationship between M1–M2 switches and final NP %
% Starting NP count
start_NP = 4009;

% Final NP counts for trials 4:33
final_NP = untreatedcontrol_NP(4:33, 1001);

% Compute % remaining
percent_remaining = (final_NP / start_NP) * 100;

% Number of switches per trial (from previous analysis)
switches = intersection_counts;  % should be length 30

% Scatter plot
figure;
scatter(switches, percent_remaining, 60, 'filled');
xlabel('Number of Switches (M1–M2 crossings)');
ylabel('Final NP % Remaining');
title('Relationship between Dominance Switches and NP Cell Survival - Untreated Control');

hold on;

%Linear regression line
coeffs = polyfit(switches, percent_remaining, 1); % linear fit
xFit = linspace(min(switches), max(switches), 100);
yFit = polyval(coeffs, xFit);
plot(xFit, yFit, 'r-', 'LineWidth', 2);
legend('Data', 'Linear Fit', 'Location', 'best');

%Compute Pearson correlation
[R, P] = corrcoef(switches, percent_remaining);
% fprintf('Correlation coefficient: %.2f\n', R(1,2));
% fprintf('p-value: %.4f\n', P(1,2));
txt = sprintf('Pearson correlation coefficient = %.2f\np = %.4f', R(1,2), P(1,2));

text(0.98, 0.98, txt, ...
    'Units', 'normalized', ...      % relative to axes
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', ...
    'FontWeight', 'bold', ...
    'FontSize', 10);
%% here is where im adding the treatment groups
% HIGH DENSITY SEALANT
%clearing variables to make room for next .mat file
clearvars -except untreatedcontrol_M1 untreatedcontrol_M2 untreatedcontrol_NP time untreatedcontrol_M1_mean_curve untreatedcontrol_M2_mean_curve untreated_m2ratio_mean_curve

%First i need to combine with the beginning of untreated data
first_M1_portion = untreatedcontrol_M1(4:33, 1:250);
first_M2_portion = untreatedcontrol_M2(4:33, 1:250);

%load in next data
load('kyannahighdensitywound2.mat'); %this data is for the 30 trials of high density sealant

sealant_NP_portion = nSTORAGE(1:30,:); 
%i want a 30x1001 double for my combined variable that i wanna plot
sealant_M1_portion = m1STORAGE(1:30,:);
sealant_M2_portion = m2STORAGE(1:30,:);

%Combing each trial
sealant_M1_combined = [first_M1_portion(:, :), sealant_M1_portion(:, :)];
sealant_M2_combined = [first_M2_portion(:, :), sealant_M2_portion(:, :)];

%Sealant 30 Trials Plotting M1 vs Time
figure;
tl = tiledlayout(1,2); % 1 row, 2 columns
nexttile(tl, 1)  % left tile

title("High Density Sealant M1 Population vs Time")
xlabel("Simulation Time (days)")
ylabel("Cell Count")
hold on; 
plot(time, sealant_M1_combined, "Color",[0.7 0.6 0.9]);
xlim([0 41.67])


%data analysis
sealant_M1_mean_curve = mean(sealant_M1_combined, 1); %1 means it is averaging across rows 
plot(time, sealant_M1_mean_curve, 'LineWidth', 2, 'Color', [0 0 0]);

%bootstrapping a 95% confidence interval of the mean
sealant_M1_data_to_boot = sealant_M1_combined(:, :); %same size we want to bootstrap
sealant_M1_ci = bootci(1000,@mean,sealant_M1_data_to_boot); 
plot(time, sealant_M1_ci(1,:), 'LineWidth', 0.75, 'Color', [0 0 0], 'LineStyle','--');
plot(time, sealant_M1_ci(2,:), 'LineWidth', 0.75, 'Color', [0 0 0], 'LineStyle','--');

% 95% CI upper and lower bounds (for shading purposes)
ci_lower = sealant_M1_ci(1,:);
ci_upper = sealant_M1_ci(2,:);

%shading confidence interval
fill([time fliplr(time)], ...
     [ci_upper fliplr(ci_lower)], ...
     'y', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

%Sealant 30 Trials Plotting M2 vs Time
% figure;
hold off
nexttile(tl, 2)  % left tile
hold on; 
plot(time, sealant_M2_combined, "Color",[0.7 0.6 0.9]);
xlim([0 41.67])
title("High Density Sealant M2 Population vs Time")
xlabel("Simulation Time (days)")
ylabel("Cell Count")


%data analysis
sealant_M2_mean_curve = mean(sealant_M2_combined, 1); %1 means it is averaging across rows 
plot(time, sealant_M2_mean_curve, 'LineWidth', 2, 'Color', [0 0 0]);

%bootstrapping a 95% confidence interval of the mean
sealant_M2_data_to_boot = sealant_M2_combined(:, :); %same size we want to bootstrap
sealant_M2_ci = bootci(1000,@mean,sealant_M2_data_to_boot); 
plot(time, sealant_M2_ci(1,:), 'LineWidth', 0.75, 'Color', [0 0 0], 'LineStyle','--');
plot(time, sealant_M2_ci(2,:), 'LineWidth', 0.75, 'Color', [0 0 0], 'LineStyle','--');

% 95% CI upper and lower bounds (for shading purposes)
ci_lower = sealant_M2_ci(1,:);
ci_upper = sealant_M2_ci(2,:);

%shading confidence interval
fill([time fliplr(time)], ...
     [ci_upper fliplr(ci_lower)], ...
     'y', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

%Plotting Means
figure;
plot(time, sealant_M1_mean_curve, 'LineWidth', 2, 'Color', 'r');
hold on
plot(time, sealant_M2_mean_curve, 'LineWidth', 2, 'Color', 'b');
title("Average M1 and M2 Populations Over Time - High Density Sealant")
xlabel("Simulation Time (days)")
ylabel("Cell Count")
legend('M1','M2')

% trying to plot m1 and m2 for each trial independently to see if we can see any oscilations
%and also counting number of change overs
intersection_counts = zeros(1, 30); %  zero vectorstore counts for each trial
all_intervals = []; %zero vector

%USING TILEDLAYOUT
figure;
tl = tiledlayout(6,5);
title(tl,'High Density Sealant M1 and M2 Populations vs Time');
xlabel(tl,'Simulation Time (days)');
ylabel(tl,'Cell Count');

sealant_M1_combined = [first_M1_portion(:, :), sealant_M1_portion(:, :)];


for i = 1:30
   
    nexttile(tl)
    
    plot(time, sealant_M1_combined(i, :), "Color",'r');
    hold on;
    plot(time, sealant_M2_combined(i, :), "Color",'b');

    title(['Trial', num2str(i)]);
    
    xlim([0 41.67])


    %Intersection stuff
    diff_signal = sealant_M1_combined(i,:) - sealant_M2_combined(i,:);

    % Find where the sign changes
    sign_changes = diff_signal(1:end-1) .* diff_signal(2:end) < 0;

    intersection_counts(i) = sum(sign_changes);
   
    %Trying to determine the time intervals between intersections
    % Convert indices to times
    t_cross = time(sign_changes);

    % Compute time intervals between crossings
    dt = diff(t_cross);

    % Store them
    all_intervals = [all_intervals, dt];

     %existing intersection count
    count = intersection_counts(i);

    % Add text annotation inside the plot
text(0.98, 0.98, ['Dominance shifts: ', num2str(count)], ...
    'Units', 'normalized', ...         % position relative to axes (0-1)
    'HorizontalAlignment', 'right', ...% align text to the right
    'VerticalAlignment', 'top', ...    % align text to the top
    'FontWeight', 'bold', ...
    'FontSize', 10, ...
    'Margin', 2);                      % small padding inside box

end
legend('M1', 'M2','Location','eastoutside');

% TRYING TO PLOT M2:M1 RATIO OVER TIME FOR ALL TRIALS
ratios_overtime = []; %initializing zero vector for thing i want to calculate/plot

%USING TILEDLAYOUT
figure;
tl = tiledlayout(6,5);
title(tl,'M2:M1 Ratio with High Density Sealant vs Time');
xlabel(tl,'Simulation Time (days)');
ylabel(tl,'M2:M1 Ratio');
dt = 60;
time_days = (0:1000) * dt / 1440;

for i = 1:30
  
    nexttile(tl)
    %calculate ratio
    ratios_overtime = sealant_M2_combined(i, :)./(sealant_M1_combined(i, :) + sealant_M2_combined(i, :)) ;
    plot(time_days, ratios_overtime, "Color",[0.859 0.169 0.224]);
  

    title(['Trial ', num2str(i)]);
    
    xlim([0 42])
    ylim([0.3 0.7])

    all_ratios(i,:) = ratios_overtime;
end

%making plot of mean curve across all trials
figure;
sealant_m2ratio_mean_curve = mean(all_ratios,1);
plot(time_days, sealant_m2ratio_mean_curve,'k','LineWidth',2, "Color",[0.859 0.169 0.224]);
xlabel("Simulation Time (days)")
ylabel("M2:M1 Ratio")
title("Mean M2:M1 Ratio vs Simulated Time - High Density Sealant")
xlim([0 42]);
ylim([0.48 0.62])
%% making a histogram of the number of intersections between M1 and M2 population numbers
figure;
histogram(intersection_counts);
title('Distribution of Intersection Counts Across High Desnity Sealant Trials');
xlabel('Number of Intersections');
ylabel('Frequency');

%making histogram of timing betwe osciallations
figure;
histogram(all_intervals);

title('Distribution of Time Between Intersections - High Desnity Sealant ');
xlabel('Simulation Time (days) Between Crossings');
ylabel('Frequency');
% Relationship between M1–M2 switches and final NP %
% Starting NP count
start_NP = 4009;

% Final NP counts for trials 4:33
final_NP = sealant_NP_portion(1:30, 751);

% Compute % remaining
percent_remaining = (final_NP / start_NP) * 100;

%Number of switches per trial
switches = intersection_counts;  % should be length 30

%Scatter plot
figure;
scatter(switches, percent_remaining, 60, 'filled');
xlabel('Number of Switches (M1–M2 crossings)');
ylabel('Final NP % Remaining');
title('Relationship between Dominance Switches and NP Cell Survival - High Desnity Sealant');

hold on;

%Linear regression line
coeffs = polyfit(switches, percent_remaining, 1); % linear fit
xFit = linspace(min(switches), max(switches), 100);
yFit = polyval(coeffs, xFit);
plot(xFit, yFit, 'r-', 'LineWidth', 2);
legend('Data', 'Linear Fit', 'Location', 'best');

%Compute Pearson correlation
[R, P] = corrcoef(switches, percent_remaining);
% fprintf('High Desnity Sealant Correlation coefficient: %.2f\n', R(1,2));
% fprintf('p-value: %.4f\n', P(1,2));
txt = sprintf('Pearson correlation coefficient = %.2f\np = %.4f', R(1,2), P(1,2));

%Add text to axes (normalized units so it stays in top-right)
text(0.98, 0.98, txt, ...
    'Units', 'normalized', ...      % relative to axes
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', ...
    'FontWeight', 'bold', ...
    'FontSize', 10);
%% M2 CELL THERAPY
%clearing variables to make room for next .mat file
clearvars -except untreatedcontrol_M1 untreatedcontrol_M2 untreatedcontrol_NP time untreatedcontrol_M1_mean_curve untreatedcontrol_M2_mean_curve untreated_m2ratio_mean_curve sealant_m2ratio_mean_curve sealant_M2_combined sealant_M1_combined sealant_NP_portion

%First i need to combine with the beginning of untreated data
first_M1_portion = untreatedcontrol_M1(4:33, 1:250);
first_M2_portion = untreatedcontrol_M2(4:33, 1:250);

%load in next data
load('kyannam2celltherapy.mat'); %this data is for the 30 trials of M2 Cell Therapy

celltherapy_NP_portion = nSTORAGE(1:30,:); 
%i want a 30x1001 double for my combined variable that i wanna plot
celltherapy_M1_portion = m1STORAGE(1:30,:);
celltherapy_M2_portion = m2STORAGE(1:30,:);

%Combing each trial
celltherapy_M1_combined = [first_M1_portion(:, :), celltherapy_M1_portion(:, :)];
celltherapy_M2_combined = [first_M2_portion(:, :), celltherapy_M2_portion(:, :)];

%celltherapy 30 Trials Plotting M1 vs Time
figure;
tl = tiledlayout(1,2); % 1 row, 2 columns
nexttile(tl, 1)  % left tile

title("M2 Cell Therapy M1 Population vs Time")
xlabel("Simulation Time (days)")
ylabel("Cell Count")
hold on; 
plot(time, celltherapy_M1_combined, "Color",[0.7 0.6 0.9]);
xlim([0 41.67])


%data analysis
celltherapy_M1_mean_curve = mean(celltherapy_M1_combined, 1); %1 means it is averaging across rows 
plot(time, celltherapy_M1_mean_curve, 'LineWidth', 2, 'Color', [0 0 0]);

%bootstrapping a 95% confidence interval of the mean
celltherapy_M1_data_to_boot = celltherapy_M1_combined(:, :); %same size we want to bootstrap
celltherapy_M1_ci = bootci(1000,@mean,celltherapy_M1_data_to_boot); 
plot(time, celltherapy_M1_ci(1,:), 'LineWidth', 0.75, 'Color', [0 0 0], 'LineStyle','--');
plot(time, celltherapy_M1_ci(2,:), 'LineWidth', 0.75, 'Color', [0 0 0], 'LineStyle','--');

% 95% CI upper and lower bounds (for shading purposes)
ci_lower = celltherapy_M1_ci(1,:);
ci_upper = celltherapy_M1_ci(2,:);

%shading confidence interval
fill([time fliplr(time)], ...
     [ci_upper fliplr(ci_lower)], ...
     'y', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

%celltherapy 30 Trials Plotting M2 vs Time
% figure;
hold off
nexttile(tl, 2)  % left tile
hold on; 
plot(time, celltherapy_M2_combined, "Color",[0.7 0.6 0.9]);
xlim([0 41.67])
title("M2 Cell Therapy M2 Population vs Time")
xlabel("Simulation Time (days)")
ylabel("Cell Count")


%data analysis
celltherapy_M2_mean_curve = mean(celltherapy_M2_combined, 1); %1 means it is averaging across rows 
plot(time, celltherapy_M2_mean_curve, 'LineWidth', 2, 'Color', [0 0 0]);

%bootstrapping a 95% confidence interval of the mean
celltherapy_M2_data_to_boot = celltherapy_M2_combined(:, :); %same size we want to bootstrap
celltherapy_M2_ci = bootci(1000,@mean,celltherapy_M2_data_to_boot); 
plot(time, celltherapy_M2_ci(1,:), 'LineWidth', 0.75, 'Color', [0 0 0], 'LineStyle','--');
plot(time, celltherapy_M2_ci(2,:), 'LineWidth', 0.75, 'Color', [0 0 0], 'LineStyle','--');

% 95% CI upper and lower bounds (for shading purposes)
ci_lower = celltherapy_M2_ci(1,:);
ci_upper = celltherapy_M2_ci(2,:);

%shading confidence interval
fill([time fliplr(time)], ...
     [ci_upper fliplr(ci_lower)], ...
     'y', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

%Plotting Means
figure;
plot(time, celltherapy_M1_mean_curve, 'LineWidth', 2, 'Color', 'r');
hold on
plot(time, celltherapy_M2_mean_curve, 'LineWidth', 2, 'Color', 'b');
title("Average M1 and M2 Populations Over Time - M2 Cell Therapy")
xlabel("Simulation Time (days)")
ylabel("Cell Count")
legend('M1','M2')

% trying to plot m1 and m2 for each trial independently to see if we can see any oscilations
%and also counting number of change overs
intersection_counts = zeros(1, 30); %  zero vectorstore counts for each trial
all_intervals = []; %zero vector

%USING TILEDLAYOUT
figure;
tl = tiledlayout(6,5);
title(tl,'M2 Cell Therapy M1 and M2 Populations vs Time');
xlabel(tl,'Simulation Time (days)');
ylabel(tl,'Cell Count');

celltherapy_M1_combined = [first_M1_portion(:, :), celltherapy_M1_portion(:, :)];


for i = 1:30
   
    nexttile(tl)
    
    plot(time, celltherapy_M1_combined(i, :), "Color",'r');
    hold on;
    plot(time, celltherapy_M2_combined(i, :), "Color",'b');

    title(['Trial', num2str(i)]);
    
    xlim([0 41.67])


    %Intersection stuff
    diff_signal = celltherapy_M1_combined(i,:) - celltherapy_M2_combined(i,:);

    % Find where the sign changes
    sign_changes = diff_signal(1:end-1) .* diff_signal(2:end) < 0;

    intersection_counts(i) = sum(sign_changes);
   
    %Trying to determine the time intervals between intersections
    % Convert indices to times
    t_cross = time(sign_changes);

    %Compute time intervals between crossings
    dt = diff(t_cross);

    %Store them
    all_intervals = [all_intervals, dt];

     %existing intersection count
    count = intersection_counts(i);

    % Add text annotation inside the plot
text(0.98, 0.98, ['Dominance shifts: ', num2str(count)], ...
    'Units', 'normalized', ...         % position relative to axes (0-1)
    'HorizontalAlignment', 'right', ...% align text to the right
    'VerticalAlignment', 'top', ...    % align text to the top
    'FontWeight', 'bold', ...
    'FontSize', 10, ...
    'Margin', 2);                      % small padding inside box

end
legend('M1', 'M2','Location','eastoutside');

% TRYING TO PLOT M2:M1 RATIO OVER TIME FOR ALL TRIALS
ratios_overtime = []; %initializing zero vector for thing i want to calculate/plot

%USING TILEDLAYOUT
figure;
tl = tiledlayout(6,5);
title(tl,'M2:M1 Ratio with M2 Cell Therapy vs Time');
xlabel(tl,'Simulation Time (days)');
ylabel(tl,'M2:M1 Ratio');
dt = 60;
time_days = (0:1000) * dt / 1440;

for i = 1:30
  
    nexttile(tl)
    %calculate ratio
    ratios_overtime = celltherapy_M2_combined(i, :)./(celltherapy_M1_combined(i, :) + celltherapy_M2_combined(i, :)) ;
    plot(time_days, ratios_overtime, "Color",[0.000 0.553 0.835]);
  

    title(['Trial ', num2str(i)]);
    
    xlim([0 42])
    ylim([0.3 0.7])

    all_ratios(i,:) = ratios_overtime;
end

%making plot of mean curve across all trials
figure;
celltherapy_m2ratio_mean_curve = mean(all_ratios,1);
plot(time_days, celltherapy_m2ratio_mean_curve,'k','LineWidth',2, "Color",[0.000 0.553 0.835]);
xlabel("Time (Simulation days)")
ylabel("M2:M1 Ratio")
title("Mean M2:M1 Ratio vs Simulated Time - M2 Cell Therapy")
xlim([0 42]);
ylim([0.48 0.62])

%% 

%making a histogram of the number of intersections
figure;
histogram(intersection_counts);
title('Distribution of Intersection Counts Across M2 Cell Therapy Trials');
xlabel('Number of Intersections');
ylabel('Frequency');

%making histogram of timing betwe osciallations
figure;
histogram(all_intervals);

title('Distribution of Time Between Intersections - M2 Cell Therapy ');
xlabel('Simulation Time (days) Between Crossings');
ylabel('Frequency');
% Relationship between M1–M2 switches and final NP %
% Starting NP count
start_NP = 4009;

% Final NP counts for trials 4:33
final_NP = celltherapy_NP_portion(1:30, 751);

% Compute % remaining
percent_remaining = (final_NP / start_NP) * 100;

% Number of switches per trial (from previous analysis)
switches = intersection_counts;  % should be length 30

% Scatter plot
figure;
scatter(switches, percent_remaining, 60, 'filled');
xlabel('Number of Switches (M1–M2 crossings)');
ylabel('Final NP % Remaining');
title('Relationship between Dominance Switches and NP Cell Survival - M2 Cell Therapy');

hold on;

% Linear regression line
coeffs = polyfit(switches, percent_remaining, 1); % linear fit
xFit = linspace(min(switches), max(switches), 100);
yFit = polyval(coeffs, xFit);
plot(xFit, yFit, 'r-', 'LineWidth', 2);
legend('Data', 'Linear Fit', 'Location', 'best');

% Compute Pearson correlation
[R, P] = corrcoef(switches, percent_remaining);
% fprintf('M2 Cell Therapy Correlation coefficient: %.2f\n', R(1,2));
% fprintf('p-value: %.4f\n', P(1,2));
% Create text string
txt = sprintf('Pearson correlation coefficient = %.2f\np = %.4f', R(1,2), P(1,2));

% Add text to axes (normalized units so it stays in top-right)
text(0.98, 0.98, txt, ...
    'Units', 'normalized', ...      % relative to axes
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', ...
    'FontWeight', 'bold', ...
    'FontSize', 10);

%% Cytokine Loaded Gel
%clearing variables to make room for next .mat file
clearvars -except untreatedcontrol_M1 untreatedcontrol_M2 untreatedcontrol_NP time untreatedcontrol_M1_mean_curve untreatedcontrol_M2_mean_curve celltherapy_m2ratio_mean_curve untreated_m2ratio_mean_curve sealant_m2ratio_mean_curve sealant_M1_combined sealant_M2_combined celltherapy_M2_combined celltherapy_M1_combined sealant_NP_portion celltherapy_NP_portion

%First i need to combine with the beginning of untreated data
first_M1_portion = untreatedcontrol_M1(4:33, 1:250);
first_M2_portion = untreatedcontrol_M2(4:33, 1:250);

%load in next data
load('kyannasoftgel.mat'); %this data is for the 30 trials of cytokine gel

cytokine_NP_portion = nSTORAGE(1:30,:); %need to play with index cuz i think i have an extra 
%i want a 30x1001 double for my combined variable that i wanna plot
cytokine_M1_portion = m1STORAGE(1:30,:);
cytokine_M2_portion = m2STORAGE(1:30,:);

%Combing each trial
cytokine_M1_combined = [first_M1_portion(:, :), cytokine_M1_portion(:, :)];
cytokine_M2_combined = [first_M2_portion(:, :), cytokine_M2_portion(:, :)];

%cytokine 30 Trials Plotting M1 vs Time
figure;
tl = tiledlayout(1,2); % 1 row, 2 columns
nexttile(tl, 1)  % left tile

title("M2 Cell Therapy M1 Population vs Time")
xlabel("Simulation Time (days)")
ylabel("Cell Count")
hold on; 
plot(time, cytokine_M1_combined, "Color",[0.7 0.6 0.9]);
xlim([0 41.67])


%data analysis
cytokine_M1_mean_curve = mean(cytokine_M1_combined, 1); %1 means it is averaging across rows 
plot(time, cytokine_M1_mean_curve, 'LineWidth', 2, 'Color', [0 0 0]);

%bootstrapping a 95% confidence interval of the mean
cytokine_M1_data_to_boot = cytokine_M1_combined(:, :); %same size we want to bootstrap
cytokine_M1_ci = bootci(1000,@mean,cytokine_M1_data_to_boot); 
plot(time, cytokine_M1_ci(1,:), 'LineWidth', 0.75, 'Color', [0 0 0], 'LineStyle','--');
plot(time, cytokine_M1_ci(2,:), 'LineWidth', 0.75, 'Color', [0 0 0], 'LineStyle','--');

% 95% CI upper and lower bounds (for shading purposes)
ci_lower = cytokine_M1_ci(1,:);
ci_upper = cytokine_M1_ci(2,:);

%shading confidence interval
fill([time fliplr(time)], ...
     [ci_upper fliplr(ci_lower)], ...
     'y', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

%cytokine 30 Trials Plotting M2 vs Time
% figure;
hold off
nexttile(tl, 2)  % left tile
hold on; 
plot(time, cytokine_M2_combined, "Color",[0.7 0.6 0.9]);
xlim([0 41.67])
title("Cytokine Loaded Gel M2 Population vs Time")
xlabel("Simulation Time (days)")
ylabel("Cell Count")


%data analysis
cytokine_M2_mean_curve = mean(cytokine_M2_combined, 1); %1 means it is averaging across rows 
plot(time, cytokine_M2_mean_curve, 'LineWidth', 2, 'Color', [0 0 0]);

%bootstrapping a 95% confidence interval of the mean
cytokine_M2_data_to_boot = cytokine_M2_combined(:, :); %same size we want to bootstrap
cytokine_M2_ci = bootci(1000,@mean,cytokine_M2_data_to_boot); 
plot(time, cytokine_M2_ci(1,:), 'LineWidth', 0.75, 'Color', [0 0 0], 'LineStyle','--');
plot(time, cytokine_M2_ci(2,:), 'LineWidth', 0.75, 'Color', [0 0 0], 'LineStyle','--');

% 95% CI upper and lower bounds (for shading purposes)
ci_lower = cytokine_M2_ci(1,:);
ci_upper = cytokine_M2_ci(2,:);

%shading confidence interval
fill([time fliplr(time)],[ci_upper fliplr(ci_lower)],'y', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

%Plotting Means
figure;
plot(time, cytokine_M1_mean_curve, 'LineWidth', 2, 'Color', 'r');
hold on
plot(time, cytokine_M2_mean_curve, 'LineWidth', 2, 'Color', 'b');
title("Average M1 and M2 Populations Over Time - Cytokine Loaded Gel")
xlabel("Simulation Time (days)")
ylabel("Cell Count")
legend('M1','M2')

% trying to plot m1 and m2 for each trial independently to see if we can see any oscilations
%and also counting number of change overs
intersection_counts = zeros(1, 30); %zero vector store counts for each trial
all_intervals = []; %zero vector

%USING TILEDLAYOUT
figure;
tl = tiledlayout(6,5);
title(tl,'Cytokine Loaded Gel M1 and M2 Populations vs Time');
xlabel(tl,'Simulation Time (days)');
ylabel(tl,'Cell Count');

cytokine_M1_combined = [first_M1_portion(:, :), cytokine_M1_portion(:, :)];


for i = 1:30
   
    nexttile(tl)
    
    plot(time, cytokine_M1_combined(i, :), "Color",'r');
    hold on;
    plot(time, cytokine_M2_combined(i, :), "Color",'b');

    title(['Trial', num2str(i)]);
    
    xlim([0 41.67])


    %Intersection stuff
    diff_signal = cytokine_M1_combined(i,:) - cytokine_M2_combined(i,:);

    %Find where the sign changes
    sign_changes = diff_signal(1:end-1) .* diff_signal(2:end) < 0;

    intersection_counts(i) = sum(sign_changes);
   
    %Trying to determine the time intervals between intersections
    %Convert indices to times
    t_cross = time(sign_changes);

    %Compute time intervals between crossings
    dt = diff(t_cross);

    %Store them
    all_intervals = [all_intervals, dt];

     %existing intersection count
    count = intersection_counts(i);

text(0.98, 0.98, ['Dominance shifts: ', num2str(count)], ...
    'Units', 'normalized', ...         
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', ...    
    'FontWeight', 'bold', ...
    'FontSize', 10, ...
    'Margin', 2);                     

end
legend('M1', 'M2','Location','eastoutside');

% TRYING TO PLOT M2:M1 RATIO OVER TIME FOR ALL TRIALS
ratios_overtime = []; %initializing zero vector for thing i want to calculate/plot

%USING TILEDLAYOUT
figure;
tl = tiledlayout(6,5);
title(tl,'M2:M1 Ratio with Cytokine Loaded Gel vs Time');
xlabel(tl,'Simulation Time (days)');
ylabel(tl,'M2:M1 Ratio');
dt = 60;
time_days = (0:1000) * dt / 1440;

for i = 1:30
  
    nexttile(tl)
    %calculate ratio
    ratios_overtime = cytokine_M2_combined(i, :)./(cytokine_M1_combined(i, :) + cytokine_M2_combined(i, :)) ;
    plot(time_days, ratios_overtime, "Color",[0.000 0.553 0.835]);
  

    title(['Trial ', num2str(i)]);
    
    xlim([0 42])
    ylim([0.3 0.7])

    all_ratios(i,:) = ratios_overtime;
end

%making plot of mean curve across all trials
figure;
cytokine_m2ratio_mean_curve = mean(all_ratios,1);
plot(time_days, cytokine_m2ratio_mean_curve,'k','LineWidth',2, "Color",[0.000 0.553 0.835]);
xlabel("Simulation Time (days)")
ylabel("M2:M1 Ratio")
title("Mean M2:M1 Ratio vs Simulated Time - Cytokine Loaded Gel")
xlim([0 42]);
ylim([0.48 0.62])

%% Figure 6a and 6b
%MAKING PLOT OF ALL MEAN M2:M1 RATIO TRAJECTORIES (SUPPORTS UPWARD ENDPOINT SHIFT WITH TREATMENT)
figure;
plot(time_days, untreated_m2ratio_mean_curve,'k','LineWidth',2,"Color",[0.231 0.208 0.380]);
hold on
plot(time_days(251:end), sealant_m2ratio_mean_curve(1, 251:1001 ),'k','LineWidth',2,"Color",[0.859 0.169 0.224]);
hold on
plot(time_days(250:end), celltherapy_m2ratio_mean_curve(1, 250:1001),'k','LineWidth',2,"Color",[0.000 0.553 0.835]);
hold on
plot(time_days(251:end), cytokine_m2ratio_mean_curve(1, 251:1001),'k','LineWidth',2,"Color",[0.475 0.710 0.600]);
xlabel("Simulation Time (days)", 'FontSize', 14);
ylabel("M2:M1 Ratio", 'FontSize', 14);
%title("Mean M2:M1 Ratio vs Simulated Time");
xlim([0 41.67])
legend("Untreated","High Density Sealant","M2 Cell Therapy","Cytokine Loaded Gel");


%MAKING PLOT FOR ENDPOINT OF M2:M1

untreated_end_ratio = untreatedcontrol_M2(4:33,1001) ./(untreatedcontrol_M1(4:33,1001) + untreatedcontrol_M2(4:33,1001));

sealant_end_ratio = sealant_M2_combined(1:30,1001) ./(sealant_M1_combined(1:30,1001) + sealant_M2_combined(1:30,1001));

cell_end_ratio = celltherapy_M2_combined(1:30,1001) ./ (celltherapy_M1_combined(1:30,1001) + celltherapy_M2_combined(1:30,1001));

cytokine_end_ratio = cytokine_M2_combined(1:30,1001) ./ (cytokine_M1_combined(1:30,1001) + cytokine_M2_combined(1:30,1001));

figure;
hold on

%organize data
groups = {untreated_end_ratio, sealant_end_ratio, cell_end_ratio, cytokine_end_ratio};

group_names = ["Untreated","High Density Sealant","M2 Cell Therapy","Cytokine Loaded Gel"];

%my color palette
colors = [
    0.231 0.208 0.380;   % untreated (dark purple)
    0.859 0.169 0.224;   % sealant
    0.000 0.553 0.835;   % cell therapy
    0.475 0.710 0.600    % cytokine
];

jitter_strength = 0.08;

for g = 1:4
    
    y = groups{g};
    
    %jitter x positions so points don’t overlap
    x = g + (rand(size(y)) - 0.5) * jitter_strength;
    
    %scatter points
    scatter(x, y,60, colors(g,:), 'filled', ...
        'MarkerFaceAlpha', 0.5)
    
    %mean
    mu = mean(y);
    
    %95% CI via bootstrapping
    ci = bootci(1000, @mean, y);
    
    %CI line
    plot([g g], [ci(1) ci(2)], 'k', 'LineWidth', 1.5)
    
    %mean point
    plot(g, mu, 'k.', 'MarkerSize', 25)
end

%formatting
xticks(1:4)
xticklabels(group_names)

ylabel('Endpoint M2:M1 Ratio', 'FontSize', 14)
%title('Endpoint M2:M1 Ratio Across Treatment Groups')

ylim([0.2 1])
xlim([0.5 4.5])
box on

%Stats for that graph

%organize data
values = [untreated_end_ratio;
          sealant_end_ratio;
          cell_end_ratio;
          cytokine_end_ratio];

group = [ones(length(untreated_end_ratio),1);
         2*ones(length(sealant_end_ratio),1);
         3*ones(length(cell_end_ratio),1);
         4*ones(length(cytokine_end_ratio),1)];

group_names = ["Untreated","High Density Sealant","M2 Cell Therapy","Cytokine Loaded Gel"];

%one-way ANOVA
[p,~,stats] = anova1(values, group, 'off');

fprintf('ANOVA p-value: %.4e\n', p);

%post-hoc comparison (Tukey)
results = multcompare(stats, 'Display','off');

%add stars comparing to control (group 1)
y_max = max(values);
y_step = 0.02;   % vertical spacing
line_count = 0;

for i = 1:size(results,1)
    
    g1 = results(i,1);
    g2 = results(i,2);
    p_val = results(i,6);
    
    % only compare to control (group 1)
    if ~(g1 == 1 || g2 == 1)
        continue
    end
    
    % significance levels
    if p_val < 0.001
        stars = '***';
    elseif p_val < 0.01
        stars = '**';
    elseif p_val < 0.05
        stars = '*';
    else
        continue
    end
    
    % determine comparison group
    other = g2;
    if g2 == 1
        other = g1;
    end
    
    line_count = line_count + 1;
    y = y_max + line_count * y_step;
    
    % draw line
    plot([1 other], [y y], 'k', 'LineWidth',1.2)
    
    % draw asterisks
    text(mean([1 other]), y + 0.005, stars, ...
        'HorizontalAlignment','center', ...
        'FontSize',18, ...
        'FontWeight','bold')
end

comparisons = [
    2 3;  
    3 4;  
    2 4  
];

y_max = max(values);
y_step = 0.04;   %increase if overlapping
offset = 0.005;  %space between line and stars

for c = 1:size(comparisons,1)
    
    g1 = comparisons(c,1);
    g2 = comparisons(c,2);
    
    %find matching row in multcompare results
    for i = 1:size(results,1)
        r1 = results(i,1);
        r2 = results(i,2);
        
        if (r1 == g1 && r2 == g2) || (r1 == g2 && r2 == g1)
            p_val = results(i,6);
            break
        end
    end
    
    %assign stars
    if p_val < 0.001
        stars = '***';
    elseif p_val < 0.01
        stars = '**';
    elseif p_val < 0.05
        stars = '*';
    else
        stars = 'ns';
    end
    
    %vertical placement of stars and lines
   if strcmp(stars, 'ns')
    y = y_max + (c + 1.2) * y_step;  %move ns down a bit
else
    y = y_max + c * y_step;
end
    %draw line
    plot([g1 g2], [y y], 'k', 'LineWidth',1.2)
    %draw stars
   if strcmp(stars, 'ns')
    text(mean([g1 g2]), y + offset + 0.015, stars, ...
        'HorizontalAlignment','center', ...
        'FontSize',14)
else
    text(mean([g1 g2]), y + offset, stars, ...
        'HorizontalAlignment','center', ...
        'FontSize',18, ...
        'FontWeight','bold')
end
    
end

%% 

untreated_ratio = untreatedcontrol_M2(4:33,:) ./(untreatedcontrol_M1(4:33,:) + untreatedcontrol_M2(4:33,:));

sealant_ratio = sealant_M2_combined(1:30,:) ./(sealant_M1_combined(1:30,:) + sealant_M2_combined(1:30,:));

celltherapy_ratio = celltherapy_M2_combined(1:30,:) ./ (celltherapy_M1_combined(1:30,:) + celltherapy_M2_combined(1:30,:));

cytokine_ratio = cytokine_M2_combined(1:30,:) ./ (cytokine_M1_combined(1:30,:) + cytokine_M2_combined(1:30,:));
     %AUC for ratio
        aucratio_untreated = trapz(untreated_ratio, 2);
        aucratio_sealant = trapz(sealant_ratio, 2);
        aucratio_cell = trapz(celltherapy_ratio, 2);
        aucratio_cytokine = trapz(cytokine_ratio, 2);
        % combining
        auc_all = [aucratio_untreated; aucratio_sealant; aucratio_cell; aucratio_cytokine];
        endpoint_np_all  = [untreatedcontrol_NP(4:33, 1001)./4009*100; sealant_NP_portion(1:30, 751)./4009*100;  celltherapy_NP_portion(1:30, 751)./4009*100;  cytokine_NP_portion(1:30, 751)./4009*100];
        %   scatter(auc_all, np_all)
       
        colors = [
        0.231 0.208 0.380;   % untreated (dark purple)
        0.859 0.169 0.224;   % sealant
        0.000 0.553 0.835;   % cell therapy
        0.475 0.710 0.600]; % cytokine
        % scatter points
        figure;
        scatter(aucratio_untreated, untreatedcontrol_NP(4:33, 1001)./4009*100, 50, [59 53 97]/255, 'filled');
        hold on
    scatter(aucratio_sealant, sealant_NP_portion(1:30, 751)./4009*100, 50, [219 43 57]/255, 'filled');
    hold on
    scatter(aucratio_cell, celltherapy_NP_portion(1:30, 751)./4009*100, 50, [0 141 213]/255, 'filled');
    hold on
    scatter(aucratio_cytokine, cytokine_NP_portion(1:30, 751)./4009*100, 50, [121 181 153]/255, 'filled');
    
        %scatter(auc_all,  endpoint_np_all, 30, colors(g,:), 'filled', 'MarkerFaceAlpha', 0.5)
         xlabel('M2:M1 Ratio AUC', 'FontSize', 14')
        ylabel('Endpoint NP % Remaining','FontSize', 14)
      %plotting line
        mdl = fitlm(auc_all, endpoint_np_all);
    
    hold on
    x = linspace(min(auc_all), max(auc_all), 100)';
    y = predict(mdl, x);
    
    plot(x, y, 'k', 'LineWidth', 2)
        %correlation
        [r, p] = corr(auc_all, endpoint_np_all);
        % Create text string
        txt = sprintf('Pearson correlation coefficient = %.2f\np = %.4f', r, p);
    
        % Add text to axes (normalized units so it stays in top-right)
        text(0.98, 0.98, txt, ...
        'Units', 'normalized', ...      % relative to axes
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'top', ...
        'FontWeight', 'bold', ...
        'FontSize', 14);
        leg = legend("Untreated","High Density Sealant","M2 Cell Therapy","Cytokine Loaded Gel");
        leg.FontSize = 12;

%DOING IT FOR INDIVIDUAL GROUPS
need_title = figure;
title("M2:M1 Ratio AUC vs Endpoint % NP Remaining by Treatment")
tiles = tiledlayout(2,2);
colors = [
        0.231 0.208 0.380;   % untreated (dark purple)
        0.859 0.169 0.224;   % sealant
        0.000 0.553 0.835;   % cell therapy
        0.475 0.710 0.600]; % cytokine
%UNTREATED
ax_new = nexttile(tiles);

scatter(aucratio_untreated, untreatedcontrol_NP(4:33,1001)./4009*100, 50, colors(1,:), 'filled');
hold on;

mdl = fitlm(aucratio_untreated, untreatedcontrol_NP(4:33,1001)./4009*100);
xfit = linspace(min(aucratio_untreated), max(aucratio_untreated), 100)';
yfit = predict(mdl, xfit);
plot(xfit, yfit, 'k', 'LineWidth', 1.5);

[r,p] = corr(aucratio_untreated, untreatedcontrol_NP(4:33,1001)./4009*100);

title(sprintf('Untreated (r=%.2f)', r),'FontSize', 14);
xlabel('M2:M1 Ratio AUC','FontSize', 14);
ylabel('Endpoint % NP Remaining','FontSize', 14);


%sealant
ax_new = nexttile(tiles);

x = aucratio_sealant;
ydata = sealant_NP_portion(1:30,751)./4009*100;

scatter(x, ydata, 50, colors(2,:), 'filled');
hold on;

mdl = fitlm(x, ydata);
xfit = linspace(min(x), max(x), 100)';
yfit = predict(mdl, xfit);
plot(xfit, yfit, 'k', 'LineWidth', 1.5);

[r,p] = corr(x, ydata);

title(sprintf('High Density Sealant (r=%.2f)', r),'FontSize', 14);
xlabel('M2:M1 Ratio AUC','FontSize', 14);
ylabel('Endpoint % NP Remaining', 'FontSize', 14);

%cell therapy
ax_new = nexttile(tiles);


x = aucratio_cell;
ydata = celltherapy_NP_portion(1:30,751)./4009*100;

scatter(x, ydata, 50, colors(3,:), 'filled');
hold on;

mdl = fitlm(x, ydata);
xfit = linspace(min(x), max(x), 100)';
yfit = predict(mdl, xfit);
plot(xfit, yfit, 'k', 'LineWidth', 1.5);

[r,p] = corr(x, ydata);

title(sprintf('M2 Cell Therapy (r=%.2f)', r),'FontSize', 14);
xlabel('M2:M1 Ratio AUC', 'FontSize', 14);
ylabel('Endpoint % NP Remaining', 'FontSize', 14);

%cytokine loaded gel
ax_new = nexttile(tiles);


x = aucratio_cytokine;
ydata = cytokine_NP_portion(1:30,751)./4009*100;

scatter(x, ydata, 50, colors(4,:), 'filled');
hold on;

mdl = fitlm(x, ydata);
xfit = linspace(min(x), max(x), 100)';
yfit = predict(mdl, xfit);
plot(xfit, yfit, 'k', 'LineWidth', 1.5);

[r,p] = corr(x, ydata);

title(sprintf('Cytokine Loaded Gel (r=%.2f)', r),'FontSize', 14);
xlabel('M2:M1 Ratio AUC','FontSize', 14);
ylabel('Endpoint % NP Remaining', 'FontSize', 14);

%% MAKING PLOTS FOR ALL MEAN TOTAL MAC COUNT TRAJECTORIES (SUPPORTS DOWNWARD ENDPOINT SHIFT WITH TREATMENT
untreated_totalmac = untreatedcontrol_M1(4:33, :) + untreatedcontrol_M2(4:33, :);
mean_untreated_totalmac = mean(untreated_totalmac);

sealant_totalmac = sealant_M1_combined + sealant_M2_combined;
mean_sealant_totalmac = mean(sealant_totalmac);

celltherapy_totalmac = celltherapy_M1_combined + celltherapy_M2_combined;
mean_celltherapy_totalmac = mean(celltherapy_totalmac);


cytokine_totalmac = cytokine_M1_combined + cytokine_M2_combined;
mean_cytokine_totalmac = mean(cytokine_totalmac);

%making AUC plot with correlations
%for each trial and for each treatment, I want to compute AUC and store
%endpoint np cell survival, then plot x = AUC  y = NP survival 

    %AUC for TOTAL MACROPHAGE
    auc_untreated = trapz(untreated_totalmac, 2);
    auc_sealant = trapz(sealant_totalmac, 2);
    auc_cell = trapz(celltherapy_totalmac, 2);
    auc_cytokine = trapz(cytokine_totalmac, 2);
    % combining
    auc_all = [auc_untreated; auc_sealant; auc_cell; auc_cytokine];
    endpoint_np_all  = [untreatedcontrol_NP(4:33, 1001)./4009*100; sealant_NP_portion(1:30, 751)./4009*100;  celltherapy_NP_portion(1:30, 751)./4009*100;  cytokine_NP_portion(1:30, 751)./4009*100];
    %   scatter(auc_all, np_all)
   
    colors = [
    0.231 0.208 0.380;   % untreated (dark purple)
    0.859 0.169 0.224;   % sealant
    0.000 0.553 0.835;   % cell therapy
    0.475 0.710 0.600]; % cytokine
    % scatter points
    figure;
    scatter(auc_untreated, untreatedcontrol_NP(4:33, 1001)./4009*100, 50, [59 53 97]/255, 'filled');
    hold on
scatter(auc_sealant, sealant_NP_portion(1:30, 751)./4009*100, 50, [219 43 57]/255, 'filled');
hold on
scatter(auc_cell, celltherapy_NP_portion(1:30, 751)./4009*100, 50, [0 141 213]/255, 'filled');
hold on
scatter(auc_cytokine, cytokine_NP_portion(1:30, 751)./4009*100, 50, [121 181 153]/255, 'filled');

    %scatter(auc_all,  endpoint_np_all, 30, colors(g,:), 'filled', 'MarkerFaceAlpha', 0.5)
     xlabel('Total Macrophage Count AUC', 'FontSize',14)
    ylabel('Endpoint NP % Remaining','FontSize', 14')
  %plotting line
    mdl = fitlm(auc_all, endpoint_np_all);

hold on
x = linspace(min(auc_all), max(auc_all), 100)';
y = predict(mdl, x);

plot(x, y, 'k', 'LineWidth', 2)
    %correlation
    [r, p] = corr(auc_all, endpoint_np_all);
    % Create text string
    txt = sprintf('Pearson correlation coefficient = %.2f\np = %.4f', r, p);

    % Add text to axes (normalized units so it stays in top-right)
    text(0.98, 0.98, txt, ...
    'Units', 'normalized', ...      % relative to axes
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', ...
    'FontWeight', 'bold', ...
    'FontSize', 12);

%TESTING IT WITHOUT CELLULAR THERPY
    %AUC
    auc_untreated = trapz(untreated_totalmac, 2);
    auc_sealant = trapz(sealant_totalmac, 2);
  
    auc_cytokine = trapz(cytokine_totalmac, 2);
    % combining
    auc_all = [auc_untreated; auc_sealant; auc_cytokine];
    endpoint_np_all  = [untreatedcontrol_NP(4:33, 1001)./4009*100; sealant_NP_portion(1:30, 751)./4009*100; cytokine_NP_portion(1:30, 751)./4009*100];
    %   scatter(auc_all, np_all)
   
    colors = [
    0.231 0.208 0.380;   % untreated (dark purple)
    0.859 0.169 0.224;   % sealant
    0.475 0.710 0.600]; % cytokine
    % scatter points
    figure;
    scatter(auc_untreated, untreatedcontrol_NP(4:33, 1001)./4009*100, 50, [59 53 97]/255, 'filled');
    hold on
scatter(auc_sealant, sealant_NP_portion(1:30, 751)./4009*100, 50, [219 43 57]/255, 'filled');
hold on
scatter(auc_cytokine, cytokine_NP_portion(1:30, 751)./4009*100, 50, [121 181 153]/255, 'filled');

    %scatter(auc_all,  endpoint_np_all, 30, colors(g,:), 'filled', 'MarkerFaceAlpha', 0.5)
     xlabel('Total Macrophage Count AUC', 'FontSize', 14)
    ylabel('Endpoint NP % Remaining','FontSize', 14)
  %plotting line
    mdl = fitlm(auc_all, endpoint_np_all);

hold on
x = linspace(min(auc_all), max(auc_all), 100)';
y = predict(mdl, x);

plot(x, y, 'k', 'LineWidth', 2)
    %correlation
    [r, p] = corr(auc_all, endpoint_np_all);
    % Create text string
    txt = sprintf('Pearson correlation coefficient = %.2f\np = %.4f', r, p);

    % Add text to axes (normalized units so it stays in top-right)
    text(0.98, 0.98, txt, ...
    'Units', 'normalized', ...      % relative to axes
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', ...
    'FontWeight', 'bold', ...
    'FontSize', 12);


%doing it dor individual treatment groups
need_title = figure;
title("Total Macrophage Count AUC vs NP Survival by Treatment")
tiles = tiledlayout(2,2);
colors = [
        0.231 0.208 0.380;   % untreated (dark purple)
        0.859 0.169 0.224;   % sealant
        0.000 0.553 0.835;   % cell therapy
        0.475 0.710 0.600]; % cytokine
%UNTREATED
ax_new = nexttile(tiles);

scatter(auc_untreated, untreatedcontrol_NP(4:33,1001)./4009*100, 50, colors(1,:), 'filled');
hold on;

mdl = fitlm(auc_untreated, untreatedcontrol_NP(4:33,1001)./4009*100);
xfit = linspace(min(auc_untreated), max(auc_untreated), 100)';
yfit = predict(mdl, xfit);
plot(xfit, yfit, 'k', 'LineWidth', 1.5);

[r,p] = corr(auc_untreated, untreatedcontrol_NP(4:33,1001)./4009*100);

title(sprintf('Untreated (r=%.2f)', r), 'FontSize',14);
xlabel('Total Macrophage Count AUC', 'FontSize', 14);
ylabel('Endpoint NP % Remaining','FontSize', 14);


%sealant
ax_new = nexttile(tiles);

x = auc_sealant;
ydata = sealant_NP_portion(1:30,751)./4009*100;

scatter(x, ydata, 50, colors(2,:), 'filled');
hold on;

mdl = fitlm(x, ydata);
xfit = linspace(min(x), max(x), 100)';
yfit = predict(mdl, xfit);
plot(xfit, yfit, 'k', 'LineWidth', 1.5);

[r,p] = corr(x, ydata);

title(sprintf('High Density Sealant (r=%.2f)', r), 'FontSize',14);
xlabel('Total Macrophage Count AUC', 'FontSize',14);
ylabel('Endpoint NP % Remaining','FontSize', 14);

%cell therapy
ax_new = nexttile(tiles);


x = auc_cell;
ydata = celltherapy_NP_portion(1:30,751)./4009*100;

scatter(x, ydata, 50, colors(3,:), 'filled');
hold on;

mdl = fitlm(x, ydata);
xfit = linspace(min(x), max(x), 100)';
yfit = predict(mdl, xfit);
plot(xfit, yfit, 'k', 'LineWidth', 1.5);

[r,p] = corr(x, ydata);

title(sprintf('M2 Cell Therapy (r=%.2f)', r), 'FontSize',14);
xlabel('Total Macrophage Count AUC','FontSize', 14);
ylabel('Endpoint NP % Remaining','FontSize', 14);

%cytokine loaded gel
ax_new = nexttile(tiles);


x = auc_cytokine;
ydata = cytokine_NP_portion(1:30,751)./4009*100;

scatter(x, ydata, 50, colors(4,:), 'filled');
hold on;

mdl = fitlm(x, ydata);
xfit = linspace(min(x), max(x), 100)';
yfit = predict(mdl, xfit);
plot(xfit, yfit, 'k', 'LineWidth', 1.5);

[r,p] = corr(x, ydata);

title(sprintf('Cytokine Loaded Gel (r=%.2f)', r), 'FontSize',14);
xlabel('Total Macrophage Count AUC','FontSize', 14);
ylabel('Endpoint NP % Remaining','FontSize', 14);

%% FIGURE 7 a & b


 %Total Macrophage Count vs Simulated Time
figure;
plot(time_days, mean_untreated_totalmac,'k','LineWidth',2,"Color",[0.231 0.208 0.380]);
hold on
plot(time_days(251:end), mean_sealant_totalmac(1, 251:1001 ),'k','LineWidth',2,"Color",[0.859 0.169 0.224]);
hold on
plot(time_days(250:end), mean_celltherapy_totalmac(1, 250:1001),'k','LineWidth',2,"Color",[0.000 0.553 0.835]);
hold on
plot(time_days(251:end), mean_cytokine_totalmac(1, 251:1001),'k','LineWidth',2,"Color",[0.475 0.710 0.600]);
xlabel("Simulation Time (days)", 'FontSize', 14);
ylabel("Mean Total Macrophage Count", 'FontSize', 14);
%title("Total Macrophage Count vs Simulated Time");

leg = legend("Untreated","High Density Sealant","M2 Cell Therapy","Cytokine Loaded Gel");
leg.FontSize = 12;

xlim([0 41.67])
%MAKING PLOT FOR ENDPOINT OF mac count FIGURE 7 b

untreated_end_maccount = (untreatedcontrol_M1(4:33,1001) + untreatedcontrol_M2(4:33,1001));

sealant_end_maccount = (sealant_M1_combined(1:30,1001) + sealant_M2_combined(1:30,1001));

cell_end_maccount = (celltherapy_M1_combined(1:30,1001) + celltherapy_M2_combined(1:30,1001));

cytokine_end_maccount =  (cytokine_M1_combined(1:30,1001) + cytokine_M2_combined(1:30,1001));

figure;
hold on

%organize data
groups = {untreated_end_maccount, sealant_end_maccount, cell_end_maccount, cytokine_end_maccount};

group_names = ["Untreated","High Density Sealant","M2 Cell Therapy","Cytokine Loaded Gel"];

%my color palette
colors = [
    0.231 0.208 0.380;   % untreated (dark purple)
    0.859 0.169 0.224;   % sealant
    0.000 0.553 0.835;   % cell therapy
    0.475 0.710 0.600    % cytokine
];

jitter_strength = 0.08;

for g = 1:4
    
    y = groups{g};
    
    %jitter x positions so points don’t overlap
    x = g + (rand(size(y)) - 0.5) * jitter_strength;
    
    %scatter points
    scatter(x, y, 65, colors(g,:), 'filled', ...
        'MarkerFaceAlpha', 0.5)
    
    %mean
    mu = mean(y);
    
    %95% CI via bootstrapping
    ci = bootci(1000, @mean, y);
    
    %CI line
    plot([g g], [ci(1) ci(2)], 'k', 'LineWidth', 1.5)
    
    %mean point
    plot(g, mu, 'k.', 'MarkerSize', 25)
end

% formatting
xticks(1:4)
xticklabels(group_names)

ylabel('Endpoint Total Macrophage Count','FontSize', 14)
%title('Endpoint Macrophage Total Across Treatment Groups')


box on

%Stats for that graph

% organize data
values = [untreated_end_maccount;
          sealant_end_maccount;
          cell_end_maccount;
          cytokine_end_maccount];

group = [ones(length(untreated_end_maccount),1);
         2*ones(length(sealant_end_maccount),1);
         3*ones(length(cell_end_maccount),1);
         4*ones(length(cytokine_end_maccount),1)];

group_names = ["Untreated","High Density Sealant","M2 Cell Therapy","Cytokine Loaded Gel"];

% one-way ANOVA
[p,~,stats] = anova1(values, group, 'off');

fprintf('ANOVA p-value: %.4e\n', p);

% post-hoc comparisons (Tukey)
results = multcompare(stats, 'Display','off');

fprintf('\nPairwise comparisons:\n');
for i = 1:size(results,1)
    g1 = results(i,1);
    g2 = results(i,2);
    p_val = results(i,6);
    
    fprintf('%s vs %s: p = %.4e\n', ...
        group_names(g1), group_names(g2), p_val);
end

%formatting the stars and lines on the figure

comparisons = [
    1 2
    1 3
    1 4
    2 3
    2 4
    3 4
];

[~,~,stats] = anova1(values, group, 'off');
results = multcompare(stats, 'Display','off');

y_max = max(values);

y_step = 0.08 * y_max;   %might need to change this number depending on data
x_offset = 0.05;         %horizontal spacing above plot

used_y = [];             %prevents collisions

for c = 1:size(comparisons,1)

    g1 = comparisons(c,1);
    g2 = comparisons(c,2);

    %find p-value
    p_val = NaN;

    for i = 1:size(results,1)
        r1 = results(i,1);
        r2 = results(i,2);

        if (r1 == g1 && r2 == g2) || (r1 == g2 && r2 == g1)
            p_val = results(i,6);
            break
        end
    end

    if isnan(p_val)
        continue
    end

    %stars
    if p_val < 0.001
        stars = '***';
    elseif p_val < 0.01
        stars = '**';
    elseif p_val < 0.05
        stars = '*';
    else
        stars = 'ns';
    end

    %assign non-overlapping y position
    y = y_max + c * y_step;

    % draw bracket
    plot([g1 g2], [y y], 'k-', 'LineWidth', 1.2)

    %lift ns slightly to avoid line cutting)
    if strcmp(stars,'ns')
        dy = 0.5*y_step;
        fontw = 'normal';
        text(mean([g1 g2]), y + dy, stars, ...
        'HorizontalAlignment','center', ...
        'FontSize',14, ...
        'FontWeight',fontw)
    else
        dy = 0.01*y_step;
        fontw = 'bold';

    text(mean([g1 g2]), y + dy, stars, ...
        'HorizontalAlignment','center', ...
        'FontSize',18, ...
        'FontWeight',fontw)
    end


end

ylim([0, y_max + (size(comparisons,1)+2)*y_step])



%% 
% 
% %making a histogram of the number of intersections
% figure;
% histogram(intersection_counts);
% title('Distribution of Intersection Counts Across Cytokine Loaded Gel Trials');
% xlabel('Number of Intersections');
% ylabel('Frequency');
% 
% %making histogram of timing betwe osciallations
% figure;
% histogram(all_intervals);
% 
% title('Distribution of Time Between Intersections - Cytokine Loaded Gel ');
% xlabel('Time Between Crossings');
% ylabel('Frequency');
% % Relationship between M1–M2 switches and final NP %
% % Starting NP count
% start_NP = 4009;
% 
% % Final NP counts for trials 4:33
% final_NP = cytokine_NP_portion(1:30, 751);
% 
% % Compute % remaining
% percent_remaining = (final_NP / start_NP) * 100;
% 
% % Number of switches per trial (from previous analysis)
% switches = intersection_counts;  % should be length 30
% 
% % Scatter plot
% figure;
% scatter(switches, percent_remaining, 60, 'filled');
% xlabel('Number of Switches (M1–M2 crossings)');
% ylabel('Final NP % Remaining');
% title('Relationship between Dominance Switches and NP Cell Survival - Cytokine Loaded Gel');
% 
% hold on;
% 
% % Linear regression line
% coeffs = polyfit(switches, percent_remaining, 1); % linear fit
% xFit = linspace(min(switches), max(switches), 100);
% yFit = polyval(coeffs, xFit);
% plot(xFit, yFit, 'r-', 'LineWidth', 2);
% legend('Data', 'Linear Fit', 'Location', 'best');
% xlim([0 11]);
% % Compute Pearson correlation
% [R, P] = corrcoef(switches, percent_remaining);
% % fprintf('Cytokine Loaded Gel Correlation coefficient: %.2f\n', R(1,2));
% % fprintf('p-value: %.4f\n', P(1,2));
% % Create text string
% txt = sprintf('Pearson correlation coefficient = %.2f\np = %.4f', R(1,2), P(1,2));
% 
% % Add text to axes (normalized units so it stays in top-right)
% text(0.98, 0.98, txt, ...
%     'Units', 'normalized', ...      % relative to axes
%     'HorizontalAlignment', 'right', ...
%     'VerticalAlignment', 'top', ...
%     'FontWeight', 'bold', ...
%     'FontSize', 10);
% 
% %trying to rectify the discreteness of it all with bucketed bins
% % --- Bin switches into groups ---
% low_idx = switches <= 2;
% mid_idx = switches >= 3 & switches <= 5;
% high_idx = switches >= 6;
% 
% % Compute means
% means = [
%     mean(percent_remaining(low_idx)), ...
%     mean(percent_remaining(mid_idx)), ...
%     mean(percent_remaining(high_idx))
% ];
% 
% % Compute standard deviations
% stds = [
%     std(percent_remaining(low_idx)), ...
%     std(percent_remaining(mid_idx)), ...
%     std(percent_remaining(high_idx))
% ];
% 
% % Plot bar chart
% figure;
% b = bar(means);
% hold on;
% 
% % Add error bars
% errorbar(1:3, means, stds, 'k.', 'LineWidth', 1.5);
% 
% % Labels
% set(gca, 'XTickLabel', {'Low (0–2)', 'Medium (3–5)', 'High (6+)'});
% xlabel('Number of Switches (Binned)');
% ylabel('Final NP % Remaining');
% title('NP Survival vs Switching (Binned)');
% 
% ylim([80 90])
% hold off;
% % Create grouping variable for ANOVA
% bin_idx = zeros(size(switches));
% 
% bin_idx(low_idx) = 1;
% bin_idx(mid_idx) = 2;
% bin_idx(high_idx) = 3;
% 
% % --- One-way ANOVA ---
% [p, tbl, stats] = anova1(percent_remaining, bin_idx, 'off');
% 
% % Display p-value
% fprintf('ANOVA p-value: %.4f\n', p);
% 
% % --- Add significance label to plot ---
% y_max = max(percent_remaining);
% 
% if p < 0.05
%     sig_label = '*';   % you could expand this to **, *** if you want
% else
%     sig_label = 'ns';  % not significant
% end
% 
% % Place text above bars
% text(2, y_max + 1, sig_label, ...
%     'HorizontalAlignment', 'center', ...
%     'FontSize', 14, ...
%     'FontWeight', 'bold');
% %text for this to write about: When trials were grouped into low, medium, and high switching regimes, a modest decrease in NP survival was observed with increasing switching. However, this trend was not statistically significant (one-way ANOVA, p = 0.1397).
% %% trying to plot all scatter plots in one figure
% % Figure numbers to group
% fig_nums = [7, 13, 19, 25];
% 
% % Create new figure with tiled layout
% figure(26);
% tl = tiledlayout(2,2);
% 
% % Store plot handles for a global legend
% all_handles = [];
% all_labels = {};
% 
% for k = 1:length(fig_nums)
%     f = figure(fig_nums(k));                 % old figure
%     axes_old = findall(f, 'type', 'axes');  % get all axes in figure
% 
%     ax_new = nexttile(tl);                   % move to next subplot
% 
%     % Copy all children and store handles for legend
%     for a = 1:length(axes_old)
%         h = copyobj(allchild(axes_old(a)), ax_new);
%         all_handles = [all_handles; h];     % store handles
%     end
% 
% 
%     % Modify subplot title
%    % Define the new titles for each subplot
%     subplot_titles = {'Untreated Control', 'High Density Sealant', 'M2 Cell Therapy', 'Cytokine Loaded Gel'};
% 
%     % Modify subplot title
%     title(ax_new, subplot_titles{k});
% 
%     close(f);  % close old figure
% end
% 
% % Global title and axis labels
% sgtitle(tl, 'Relationship between Macrophage Dominance Switches and NP Cell Survival');
% xlabel(tl, 'Number of Macrophage Dominance Switches');
% ylabel(tl, 'Endpoint NP % Remaining');

% %% =========================================================
% % MACROPHAGE DOMINANCE ANALYSIS
% % % TIME M2 DOMINANT + CORRELATION WITH NP SURVIVAL
% %first step is reload/rebuild combined matrices
% 
% load('kyannahighdensitywound2.mat');
% sealant_M1_portion = m1STORAGE(1:30,:);
% sealant_M2_portion = m2STORAGE(1:30,:);
% sealant_NP_portion = nSTORAGE(1:30,:);
% 
% load('kyannam2celltherapy.mat');
% celltherapy_M1_portion = m1STORAGE(1:30,:);
% celltherapy_M2_portion = m2STORAGE(1:30,:);
% celltherapy_NP_portion = nSTORAGE(1:30,:);
% 
% load('kyannasoftgel.mat');
% cytokine_M1_portion = m1STORAGE(1:30,:);
% cytokine_M2_portion = m2STORAGE(1:30,:);
% cytokine_NP_portion = nSTORAGE(1:30,:);
% 
% first_M1_portion = untreatedcontrol_M1(4:33,1:250);
% first_M2_portion = untreatedcontrol_M2(4:33,1:250);
% 
% sealant_M1_combined = [first_M1_portion, sealant_M1_portion];
% sealant_M2_combined = [first_M2_portion, sealant_M2_portion];
% 
% celltherapy_M1_combined = [first_M1_portion, celltherapy_M1_portion];
% celltherapy_M2_combined = [first_M2_portion, celltherapy_M2_portion];
% 
% cytokine_M1_combined = [first_M1_portion, cytokine_M1_portion];
% cytokine_M2_combined = [first_M2_portion, cytokine_M2_portion];
% 
% dt = time(2) - time(1);  % not strictly needed, but kept for clarity
% 
% % FUNCTIONAL IDEA: M2 dominant = M2 > M1 at each timepoint
% 
% % --- UNTREATED CONTROL ---
% untreated_M2_dom = untreatedcontrol_M2(4:33,:) > untreatedcontrol_M1(4:33,:);
% percent_M2_untreated = mean(untreated_M2_dom, 2) * 100;
% 
% final_NP_untreated = untreatedcontrol_NP(4:33, end);
% start_NP = untreatedcontrol_NP(4:33, 1);
% percent_remaining_untreated = (final_NP_untreated ./ start_NP) * 100;
% 
% % --- HIGH DENSITY SEALANT ---
% sealant_M2_dom = sealant_M2_combined > sealant_M1_combined;
% percent_M2_sealant = mean(sealant_M2_dom, 2) * 100;
% 
% final_NP_sealant = sealant_NP_portion(1:30, end);
% percent_remaining_sealant = (final_NP_sealant ./ start_NP(1:30)) * 100;
% 
% % --- M2 CELL THERAPY ---
% celltherapy_M2_dom = celltherapy_M2_combined > celltherapy_M1_combined;
% percent_M2_celltherapy = mean(celltherapy_M2_dom, 2) * 100;
% 
% final_NP_celltherapy = celltherapy_NP_portion(1:30, end);
% percent_remaining_celltherapy = (final_NP_celltherapy ./ start_NP(1:30)) * 100;
% 
% % --- CYTOKINE GEL ---
% cytokine_M2_dom = cytokine_M2_combined > cytokine_M1_combined;
% percent_M2_cytokine = mean(cytokine_M2_dom, 2) * 100;
% 
% final_NP_cytokine = cytokine_NP_portion(1:30, end);
% percent_remaining_cytokine = (final_NP_cytokine ./ start_NP(1:30)) * 100;
% 
% 
% all_percent_M2 = [
%     percent_M2_untreated;
%     percent_M2_sealant;
%     percent_M2_celltherapy;
%     percent_M2_cytokine
% ];
% 
% all_percent_NP = [
%     percent_remaining_untreated;
%     percent_remaining_sealant;
%     percent_remaining_celltherapy;
%     percent_remaining_cytokine
% ];
% 
% group_labels = [ ...
%     ones(size(percent_M2_untreated));
%     2*ones(size(percent_M2_sealant));
%     3*ones(size(percent_M2_celltherapy));
%     4*ones(size(percent_M2_cytokine))
% ];
% 
% % =========================
% % SCATTER PLOT: % M2 DOMINANCE vs NP SURVIVAL
% % =========================
% 
% figure;
% hold on;
% 
% scatter(percent_M2_untreated, percent_remaining_untreated, 50, [59 53 97]/255, 'filled');
% scatter(percent_M2_sealant, percent_remaining_sealant, 50, [219 43 57]/255, 'filled');
% scatter(percent_M2_celltherapy, percent_remaining_celltherapy, 50, [0 141 213]/255, 'filled');
% scatter(percent_M2_cytokine, percent_remaining_cytokine, 50, [121 181 153]/255, 'filled');
% 
% xlabel('% Time M2 Dominant');
% ylabel('Endpoint NP % Remaining');
% title('% M2 Macrophage Dominance vs NP Survival');
% 
% 
% legend({'Untreated','High Density Sealant','M2 Cell Therapy','Cytokine Loaded Gel'}, ...
%     'Location','best');
% 
% % =========================
% % CORRELATIONS (OVERALL)
% % =========================
% 
% x_all = [
%     percent_M2_untreated;
%     percent_M2_sealant;
%     percent_M2_celltherapy;
%     percent_M2_cytokine
% ];
% 
% y_all = [
%     percent_remaining_untreated;
%     percent_remaining_sealant;
%     percent_remaining_celltherapy;
%     percent_remaining_cytokine
% ];
% 
% [R, P] = corrcoef(x_all, y_all);
% 
% txt = sprintf('Overall correlation R = %.2f\np = %.4f', R(1,2), P(1,2));
% 
% text(0.98, 0.98, txt, ...
%     'Units', 'normalized', ...
%     'HorizontalAlignment', 'right', ...
%     'VerticalAlignment', 'top', ...
%     'FontWeight', 'bold');
% 
% % =========================
% % OPTIONAL: GROUPED MEANS BAR PLOT
% % =========================
% 
% figure;
% means = [
%     mean(percent_M2_untreated), ...
%     mean(percent_M2_sealant), ...
%     mean(percent_M2_celltherapy), ...
%     mean(percent_M2_cytokine)
% ];
% 
% bar(means);
% set(gca, 'XTickLabel', {'Untreated','Sealant','Cell Therapy','Cytokine'});
% ylabel('% Time M2 Dominant');
% title('Macrophage Polarization Bias Across Treatments');
% 
% % =========================
% % SCATTER PLOT: % M2:M1 Ratio vs NP SURVIVAL
% % =========================
% untreated_end_ratio = untreatedcontrol_M2(4:33,1001) ./(untreatedcontrol_M1(4:33,1001) + untreatedcontrol_M2(4:33,1001));
% 
% sealant_end_ratio = sealant_M2_combined(1:30,1001) ./(sealant_M1_combined(1:30,1001) + sealant_M2_combined(1:30,1001));
% 
% cell_end_ratio = celltherapy_M2_combined(1:30,1001) ./ (celltherapy_M1_combined(1:30,1001) + celltherapy_M2_combined(1:30,1001));
% 
% cytokine_end_ratio = cytokine_M2_combined(1:30,1001) ./ (cytokine_M1_combined(1:30,1001) + cytokine_M2_combined(1:30,1001));
% 
% 
% figure;
% hold on;
% 
% scatter(untreated_end_ratio, percent_remaining_untreated, 50, [59 53 97]/255, 'filled');
% scatter(sealant_end_ratio, percent_remaining_sealant, 50, [219 43 57]/255, 'filled');
% scatter(cell_end_ratio , percent_remaining_celltherapy, 50, [0 141 213]/255, 'filled');
% scatter(cytokine_end_ratio, percent_remaining_cytokine, 50, [121 181 153]/255, 'filled');
% 
% xlabel('% M2:M1 Ratio');
% ylabel('Endpoint NP % Remaining');
% title('% Endpoint M2:M1 Ratio  vs NP Survival');
% 
% 
% legend({'Untreated','High Density Sealant','M2 Cell Therapy','Cytokine Loaded Gel'}, ...
%     'Location','best');
% 
% % =========================
% % CORRELATIONS (OVERALL)
% % =========================
% 
% x_all = [
%     untreated_end_ratio;
%     sealant_end_ratio;
%     cell_end_ratio;
%     cytokine_end_ratio
% ];
% 
% y_all = [
%     percent_remaining_untreated;
%     percent_remaining_sealant;
%     percent_remaining_celltherapy;
%     percent_remaining_cytokine
% ];
% 
% [R, P] = corrcoef(x_all, y_all);
% 
% txt = sprintf('Overall correlation R = %.2f\np = %.4f', R(1,2), P(1,2));
% 
% text(0.98, 0.98, txt, ...
%     'Units', 'normalized', ...
%     'HorizontalAlignment', 'right', ...
%     'VerticalAlignment', 'top', ...
%     'FontWeight', 'bold');
