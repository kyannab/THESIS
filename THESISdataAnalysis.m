%this is going to be my workspace for attempting data analysis
%first goal is to plot each NP population vs time graph

%% UNTREATED CONTROL TRIALS
clear all
load('kyanna30trialsdata.mat'); %this data is for the 30 trials data, representing the untreated controls

untreatedcontrol_NP = nSTORAGE(:,:);
untreatedcontrol_M1 = m1STORAGE(:,:);
untreatedcontrol_M2 = m2STORAGE(:,:);
%% 
close all
%for each storage variable, I need to plot starting from row 4->row 33 (for untreated control because of how the saving went)

%Plotting NP population vs time
%time = linspace(0, 1001, 1001); %this is if I want time to be in terms of timesteps (data saved every 60 seconds of the 60,000 minute simulation, makes for 1000 time steps
time = (0:1000) * 60 / 1440; %this converts time to simulated 41.67 days for 60,000 minutes
figure;

title("Untreated Control Nucleus Pulposus Population vs Time")
xlabel("Simulation Time (days)")
ylabel("Cell Count")
hold on; 
for i = 4:33
    plot(time, untreatedcontrol_NP(i, :), "Color",[0.231 0.208 0.380]);
end
xlim([0 41.67])
ylim([3000 4200])


%data analysis
%bootstrapping w/ 95% CI 
untreatedcontrol_NP_mean_curve = mean(untreatedcontrol_NP, 1); %1 means it is averaging across rows 
plot(time, untreatedcontrol_NP_mean_curve, 'LineWidth', 2, 'Color', [0 0 0]);

%bootstrapping a 95% confidence interval of the mean
untreatedcontrol_NP_data_to_boot = untreatedcontrol_NP(4:33, :); %same size we want to bootstrap
untreatedcontrol_NP_ci = bootci(1000,@mean,untreatedcontrol_NP_data_to_boot); 
plot(time, untreatedcontrol_NP_ci(1,:), 'LineWidth', 0.75, 'Color', [0 0 0], 'LineStyle','--');
plot(time, untreatedcontrol_NP_ci(2,:), 'LineWidth', 0.75, 'Color', [0 0 0], 'LineStyle','--');

% 95% CI upper and lower bounds (for shading purposes)
ci_lower = untreatedcontrol_NP_ci(1,:);
ci_upper = untreatedcontrol_NP_ci(2,:);

%shading confidence interval
fill([time fliplr(time)], ...
     [ci_upper fliplr(ci_lower)], ...
    [234, 82, 111]./255, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% exportgraphics(gcf,'testresolution.png','Resolution',300) %seeing if export is clearer -> it is!


%% HIGH DENSITY SEALANT
%clearing variables to make room for next .mat file
clearvars -except untreatedcontrol_M1 untreatedcontrol_M2 untreatedcontrol_NP time untreatedcontrol_NP_mean_curve

%First i need to combine with the beginning of untreated data
first_NP_portion = untreatedcontrol_NP(4:33, 1:250);

%load in next data
load('kyannahighdensitywound2.mat'); %this data is for the 30 trials of high density sealant

sealant_NP_portion = nSTORAGE(1:30,:); %need to specify 1:30 because I have an extra 
sealant_M1_portion = m1STORAGE(1:30,:);
sealant_M2_portion = m2STORAGE(1:30,:);

%plotting each trial
sealant_NP_combined = [first_NP_portion(:, :), sealant_NP_portion(:, :)];

figure;
title("High Density Sealant Nucleus Pulposus Population vs Time")
xlabel("Simulation Time (days)")
ylabel("Cell Count")
hold on; 
plot(time, sealant_NP_combined, "Color",[0.231 0.208 0.380]);
xlim([0 41.67])
ylim([3000 4200])

%Data analysis
sealant_NP_mean_curve = mean(sealant_NP_combined, 1); %1 means it is averaging across rows 
plot(time, sealant_NP_mean_curve, 'LineWidth', 2, 'Color', [0 0 0]);

%Bootstrapping a 95% confidence interval of the mean
sealant_NP_data_to_boot = sealant_NP_combined(:, :); %same size we want to bootstrap
sealant_NP_ci = bootci(1000,@mean,sealant_NP_data_to_boot); 
plot(time, sealant_NP_ci(1,:), 'LineWidth', 0.75, 'Color', [0 0 0], 'LineStyle','--');
plot(time, sealant_NP_ci(2,:), 'LineWidth', 0.75, 'Color', [0 0 0], 'LineStyle','--');

% 95% CI upper and lower bounds (for shading purposes)
ci_lower = sealant_NP_ci(1,:);
ci_upper = sealant_NP_ci(2,:);

%shading confidence interval
fill([time fliplr(time)], ...
     [ci_upper fliplr(ci_lower)], ...
     [234, 82, 111]./255, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
%% %% M2 CELL THERAPY
%clearing variables to make room for next .mat file
clearvars -except untreatedcontrol_M1 untreatedcontrol_M2 untreatedcontrol_NP time first_NP_portion untreatedcontrol_NP_mean_curve sealant_NP_mean_curve sealant_NP_combined

%load in cell therapy data
load('kyannam2celltherapy.mat'); %this data is for the 30 trials of m2 cell therapy
celltherapy_NP_portion = nSTORAGE(1:30,:); 
%i want a 30x1001 double for my combined variable that i wanna plot
celltherapy_M1_portion = m1STORAGE(1:30,:);
celltherapy_M2_portion = m2STORAGE(1:30,:);

%plotting each trial
celltherapy_NP_combined = [first_NP_portion(:, :), celltherapy_NP_portion(:, :)];

figure;
title("M2 Cell Therapy Nucleus Pulposus Population vs Time")
xlabel("Simulation Time (days)")
ylabel("Cell Count")
hold on; 
plot(time, celltherapy_NP_combined, "Color",[0.231 0.208 0.380]);
xlim([0 41.67])
ylim([3000 4200])

%data analysis
celltherapy_NP_mean_curve = mean(celltherapy_NP_combined, 1); %1 means it is averaging across rows 
plot(time, celltherapy_NP_mean_curve, 'LineWidth', 2, 'Color', [0 0 0]);

%bootstrapping a 95% confidence interval of the mean
celltherapy_NP_data_to_boot = celltherapy_NP_combined(:, :); %same size we want to bootstrap
celltherapy_NP_ci = bootci(1000,@mean,celltherapy_NP_data_to_boot); 
plot(time, celltherapy_NP_ci(1,:), 'LineWidth', 0.75, 'Color', [0 0 0], 'LineStyle','--');
plot(time, celltherapy_NP_ci(2,:), 'LineWidth', 0.75, 'Color', [0 0 0], 'LineStyle','--');

% 95% CI upper and lower bounds (for shading purposes)
ci_lower = celltherapy_NP_ci(1,:);
ci_upper = celltherapy_NP_ci(2,:);

%shading confidence interval
fill([time fliplr(time)], ...
     [ci_upper fliplr(ci_lower)], ...
     [234, 82, 111]./255, 'FaceAlpha', 0.3, 'EdgeColor', 'none');

%% %% %% Cytokine Loaded Gel
%clearing variables to make room for next .mat file
clearvars -except time untreatedcontrol_M1 untreatedcontrol_M2 untreatedcontrol_NP time first_NP_portion untreatedcontrol_NP_mean_curve sealant_NP_mean_curve celltherapy_NP_mean_curve sealant_NP_combined celltherapy_NP_combined

%load in cytokine loaded gel data
load('kyannasoftgel.mat'); %this data is for the 30 trials of cytokine loaded gel
cytokine_NP_portion = nSTORAGE(1:30,:); 
%i want a 30x1001 double for my combined variable that i wanna plot
cytokine_M1_portion = m1STORAGE(1:30,:);
cytokine_M2_portion = m2STORAGE(1:30,:);

%plotting each trial
cytokine_NP_combined = [first_NP_portion(:, :), cytokine_NP_portion(:, :)];

figure;
title("Cytokine Loaded Gel Nucleus Pulposus Population vs Time")
xlabel("Simulation Time (days)")
ylabel("Cell Count")
hold on; 
plot(time, cytokine_NP_combined, "Color",[0.231 0.208 0.380]);
xlim([0 41.67])
ylim([3000 4200])

%data analysis
cytokine_NP_mean_curve = mean(cytokine_NP_combined, 1); %1 means it is averaging across rows 
plot(time, cytokine_NP_mean_curve, 'LineWidth', 2, 'Color', [0 0 0]);

%bootstrapping a 95% confidence interval of the mean
cytokine_NP_data_to_boot = cytokine_NP_combined(:, :); %same size we want to bootstrap
cytokine_NP_ci = bootci(1000,@mean,cytokine_NP_data_to_boot); 
plot(time, cytokine_NP_ci(1,:), 'LineWidth', 0.75, 'Color', [0 0 0], 'LineStyle','--');
plot(time, cytokine_NP_ci(2,:), 'LineWidth', 0.75, 'Color', [0 0 0], 'LineStyle','--');

% 95% CI upper and lower bounds (for shading purposes)
ci_lower = cytokine_NP_ci(1,:);
ci_upper = cytokine_NP_ci(2,:);

%shading confidence interval
fill([time fliplr(time)], ...
     [ci_upper fliplr(ci_lower)], ...
     [234, 82, 111]./255, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
%% %% %% Healthy Control
%clearing variables to make room for next .mat file
clearvars -except time untreatedcontrol_M1 untreatedcontrol_M2 untreatedcontrol_NP time first_NP_portion untreatedcontrol_NP_mean_curve sealant_NP_mean_curve celltherapy_NP_mean_curve sealant_NP_combined celltherapy_NP_combined cytokine_NP_combined cytokine_NP_mean_curve
load('healthycontrol.mat');
healthycontrol_NP = nSTORAGE(:,:);

%Plotting NP population vs time
figure; 
title("Healthy Control Nucleus Pulposus Population vs Time")
xlabel("Simulation Time (days)")
ylabel("Cell Count")
hold on; 
for i = 1:30
    plot(time, healthycontrol_NP(i, :), "Color",[0.231 0.208 0.380]);
end
xlim([0 41.67])
ylim([3800 4800])

%data analysis
%i am going to try plotting average/standard deviation/standard error AND
%bootstrapping w/ 95% CI to see what the differnce is
healthycontrol_NP_mean_curve = mean(healthycontrol_NP, 1); %1 means it is averaging across rows 
plot(time, healthycontrol_NP_mean_curve, 'LineWidth', 2, 'Color', [0 0 0]);

%bootstrapping a 95% confidence interval of the mean
healthycontrol_NP_data_to_boot = healthycontrol_NP(:, :); %same size we want to bootstrap
healthycontrol_NP_ci = bootci(1000,@mean,healthycontrol_NP_data_to_boot); 
plot(time, healthycontrol_NP_ci(1,:), 'LineWidth', 0.75, 'Color', [0 0 0], 'LineStyle','--');
plot(time, healthycontrol_NP_ci(2,:), 'LineWidth', 0.75, 'Color', [0 0 0], 'LineStyle','--');

% 95% CI upper and lower bounds (for shading purposes)
ci_lower = healthycontrol_NP_ci(1,:);
ci_upper = healthycontrol_NP_ci(2,:);

%shading confidence interval
fill([time fliplr(time)], ...
     [ci_upper fliplr(ci_lower)], ...
     [234, 82, 111]./255, 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% %% TRYING TO GRAPH %degradation this was what i did all by myself
% %i want to take the last value of the combined data, average across 30 trials and then divide by 4009

%% The purpose of this section is to plot a bar graph of the endpoint NP percentage across all 30 trials for each treatment group

%final NP count across all 30 trials for each treatment group
%these values are the percentages remaining 
untreated_vals = untreatedcontrol_NP(4:33,1001) / 4009 * 100;
sealant_vals = sealant_NP_combined(:,1001) / 4009 * 100;
celltherapy_vals = celltherapy_NP_combined(:,1001) / 4009 * 100;
cytokine_vals = cytokine_NP_combined(:,1001) / 4009 * 100;
healthy_vals = healthycontrol_NP(:,1001) / 4009 * 100;

%calculating the mean and standard error of the means for %NP remaining
means = [mean(healthy_vals), mean(untreated_vals), mean(sealant_vals), mean(celltherapy_vals), mean(cytokine_vals)];

sem = [std(untreated_vals)/sqrt(30), std(healthy_vals)/sqrt(30), std(sealant_vals)/sqrt(30), std(celltherapy_vals)/sqrt(30),std(cytokine_vals)/sqrt(30)];

%creating the bar graph
figure
X = categorical({'Healthy Control','Untreated Control', 'High Density Sealant','M2 Cell Therapy', 'Cytokine Loaded Gel'});
X = reordercats(X, {'Healthy Control','Untreated Control', 'High Density Sealant','M2 Cell Therapy', 'Cytokine Loaded Gel'});

b = bar(X, means);
hold on

% Color bars
b.FaceColor = 'flat';
b.CData(1,:) = [242, 220, 93]./255; % yellow
b.CData(2,:) = [0.231 0.208 0.380];   % purple
b.CData(3,:) = [0.859 0.169 0.224];   % red
b.CData(4,:) = [0.000 0.553 0.835];   % blue
b.CData(5,:) = [0.475 0.710 0.600];   % green
ylim([80 120])

% Error bars (standard error of the mean)
errorbar(1:5, means, sem, 'k', 'linestyle','none','LineWidth',1.5)

ylabel('Endpoint %NP Remaining')
title('Average Percentage of NP Remaining Aross all Treatment Groups')
ylim([80 125])

%calculating the stats for the bar graph
all_data = [healthy_vals; untreated_vals; sealant_vals; celltherapy_vals; cytokine_vals];

group_labels = [repmat({'Healthy Control'},30,1); repmat({'Untreated Control'},30,1);
repmat({'Sealant'},30,1);
repmat({'Cell Therapy'},30,1);
repmat({'Cytokine'},30,1)];

%one-way ANOVA
[p, tbl, stats] = anova1(all_data, group_labels, 'off');

%Tukey's post hoc comparisons
results = multcompare(stats, 'Display','off');

%Adding significance stars
y_max = max(means + sem);
spacing = 1.1; % vertical spacing between lines
counter = 1;

for i = 1:size(results,1)
g1 = results(i,1);
g2 = results(i,2);
p_val = results(i,6);

if p_val < 0.05
    % Assign stars
    if p_val < 0.001
        stars = '***';
    elseif p_val < 0.01
        stars = '**';
    else
        stars = '*';
    end

    % Position
    y = y_max + counter*spacing;

    % Draw line
    plot([g1 g2], [y y], 'k', 'LineWidth', 1.5)

    % Add stars
    text(mean([g1 g2]), y + 0.2, stars, ...
        'HorizontalAlignment','center', 'FontSize', 18)

    counter = counter + 1;
end

end
%% This is my Figure 5B (Bar graph of endpoint %NP remaining, not inclusive of healthy control
%final NP count across all 30 trials for each treatment group
%these values are the percentages remaining
untreated_vals = untreatedcontrol_NP(4:33,1001) / 4009 * 100;
sealant_vals = sealant_NP_combined(:,1001) / 4009 * 100;
celltherapy_vals = celltherapy_NP_combined(:,1001) / 4009 * 100;
cytokine_vals = cytokine_NP_combined(:,1001) / 4009 * 100;

%calculating the means as above
means = [mean(untreated_vals), mean(sealant_vals), mean(celltherapy_vals), mean(cytokine_vals)];

%Bootstrap 95% CI for error bars
ci_untreated = bootci(1000,@mean,untreated_vals);
ci_sealant = bootci(1000,@mean,sealant_vals);
ci_celltherapy = bootci(1000,@mean,celltherapy_vals);
ci_cytokine = bootci(1000,@mean,cytokine_vals);

ci_lower = [ci_untreated(1), ci_sealant(1), ci_celltherapy(1), ci_cytokine(1)];

ci_upper = [ci_untreated(2), ci_sealant(2), ci_celltherapy(2), ci_cytokine(2)];

    % Convert to asymmetric error bars
    lower_err = means - ci_lower;
    upper_err = ci_upper - means;

%Making the bar graph
figure
X = categorical({'Untreated Control', 'High Density Sealant','M2 Cell Therapy', 'Cytokine Loaded Gel'});
X = reordercats(X, {'Untreated Control', 'High Density Sealant','M2 Cell Therapy', 'Cytokine Loaded Gel'});

b = bar(X, means);
hold on
ax = gca;
ax.FontSize = 12; 

% Color bars
b.FaceColor = 'flat';
b.CData(1,:) = [0.231 0.208 0.380];   % black [0.231 0.208 0.380]
b.CData(2,:) = [0.859 0.169 0.224];   % red 
b.CData(3,:) = [0.000 0.553 0.835];   % blue
b.CData(4,:) = [0.475 0.710 0.600];   % green
ylim([80 90])
% Add 95% CI error bars (this is what the current version is in the paper)
%errorbar(1:4, means, lower_err, upper_err, 'k', 'linestyle','none','LineWidth',1.5)

%Calculating standard error of mean
sem = [std(untreated_vals)/sqrt(30), std(sealant_vals)/sqrt(30), std(celltherapy_vals)/sqrt(30),std(cytokine_vals)/sqrt(30)];
% Error bars (standard error of mean)
errorbar(1:4, means, sem, 'k', 'linestyle','none','LineWidth',1.5)

ylabel('Endpoint % NP Remaining', 'FontSize', 12)
%title('Percent NP Degradation with Error Bars as 95% Confidence Intervals', 'FontSize', 18)
ylim([80 95])

%now testing stats
%one-way ANOVA (tests if any group is different overall)
all_data = [untreated_vals; sealant_vals; celltherapy_vals; cytokine_vals];

group_labels = [repmat({'Control'},30,1);
repmat({'Sealant'},30,1);
repmat({'CellTherapy'},30,1);
repmat({'Cytokine'},30,1)];

[p_anova, tbl, stats] = anova1(all_data, group_labels, 'off');

%display ANOVA result (0 means highly significant)
fprintf('One-way ANOVA p-value: %.4f\n', p_anova);

%Tukey test to compare each group
results = multcompare(stats, 'Display','off');

%adding signficance stars
hold on

y_max = max(means + upper_err); % top of CI bars
spacing = 0.8; % vertical spacing between significance lines
counter = 1;

for i = 1:size(results,1)


g1 = results(i,1);
g2 = results(i,2);
p_val = results(i,6);

% Only plot significant comparisons
if p_val < 0.05

    % Assign star level
    if p_val < 0.001
        stars = '***';
    elseif p_val < 0.01
        stars = '**';
    else
        stars = '*';
    end

    % Set vertical position
    y = y_max + counter*spacing;

    % Draw horizontal line
    plot([g1 g2], [y y], 'k', 'LineWidth', 1.5)

    % Add stars text
    text(mean([g1 g2]), y + 0.2, stars,'HorizontalAlignment','center', 'FontSize', 24)

    counter = counter + 1;
end


end
%% MOVING ON TO MACROPHAGE OPTIMIZATION
%I have data in 3 different places 
% kyannadata.mat (1 trial only, no 0.5)
% kyannarepeatdata.mat (2 trials, inclusive of 0.5)
% additionaltrials.mat (additional trial of 0.5)
%GOAL: make heatmap and do sensitivity analysis
clear all
load('kyannadata.mat');
firstfile = nSTORAGE(:,:);
trials1 = table(trial_tracker, firstfile);

clearvars -except firstfile trials1
load('kyannarepeatdata.mat')
secondfile = nSTORAGE(:,:);
trials2 = table(trial_tracker, secondfile);

clearvars -except firstfile trials1 secondfile trials2
load('additionaltrials.mat')
thirdfile = nSTORAGE(:,:);
trials3 = table(trial_tracker, thirdfile);

%the first step in organization is adding in the missing 0.5 datapoints(additionaltrials.mat)to firstfile.mat
%a summary of secondfile organization is given below; secondfile is what will be used for all processing

%To insert in the middle of the file, I need to split it. k represents the target position I want it to be in firstfile
    k = 22; % Target position where I want it to land (0.5 with 400) 
 firstfile = [firstfile(1:k-1, :); thirdfile(4, :); firstfile(k:end, :)];

    k = 16; % Target position where I want it to land (0.5 with 300) 
 firstfile = [firstfile(1:k-1, :); thirdfile(3, :); firstfile(k:end, :)];

    k = 10; % Target position where I want it to land (0.5 with 200)
 firstfile = [firstfile(1:k-1, :); thirdfile(2, :); firstfile(k:end, :)];

    k = 4; % Target position where u want it to land (0.5 with 100)
 firstfile = [firstfile(1:k-1, :); thirdfile(1, :); firstfile(k:end, :)];

%Now I am adding firstfile into secondfile (where the two repeats were stored)
%Copy last row of firstfile and add it to the end of secondfile
    secondfile(end+1, :) = firstfile(29, :); 

%now loop to get firstfile into secondfile
%starting after 0 0 is inserted at bottom indices
k = 55;  % initial insertion point in secondfile
j = 28;  % initial row in firstfile

num_insertions = 28; % total number of insertions

for iter = 1:num_insertions
    secondfile = [secondfile(1:k-1, :); firstfile(j, :); secondfile(k:end, :)];
    k = k - 2;  % decrement k by 2 each iteration
    j = j - 1;  % decrement j by 1 each iteration
end 

% summary of second file organization:
% 3 consecutive are 3 repeats for same point 
% 1:21 are all for 100 macrophages total -> 
% m2 proportion 0 (1:3), 0.2 (4:6), 0.4 (7:9), 0.5 (10:12), 0.6 (13:15), 0.8 (16:18), 1 (19:21)
% 22:42 are all 200 macrophages ->
% m2 proportion 0 (22:24), 0.2 (25:27), 0.4 (28:30), 0.5 (31:33), 0.6 (34:36),
% 0.8 (37:39), 1 (40:42)
% 43:63 are all 300 macrophages
% 64:84 are all 400 macrophages
% 85:87 is 0 m2 proportion with 0 macrophages

%NOW I AM DOING THE SENSITIVITY ANALYSIS 
%I want to see which variable impacts % Remaining Nucleus Pulposus at Endpoint the most (either total macrophage count or ratio)

%Parameters
ratios = [0, 0.2, 0.4, 0.5, 0.6, 0.8, 1];  %M2:M1 proportions (Y-axis)
totals = [100, 200, 300, 400]; %Total macrophage counts (X-axis)
num_repeats = 3; %3 repeats per parameter combo
num_ratios = length(ratios);
num_totals = length(totals);
initial_count = 4009; %starting nucleus pulposus count

%Combined data is stored in 'secondfile' (87x1001)

%Initialize matrix for % remaining
percent_remaining = nan(num_ratios, num_totals);  % rows = ratios, columnss = totals

%Loop over totals and ratios
for t_idx = 1:num_totals
    for r_idx = 1:num_ratios
        % Row indices in secondfile for this total x ratio
        start_row = (t_idx-1)*num_ratios*num_repeats + (r_idx-1)*num_repeats + 1;
        end_row = start_row + num_repeats - 1;
        
        %Extract repeat data
        repeat_data = secondfile(start_row:end_row, :);
        
        %Compute mean of final time point across repeats
        final_mean = mean(repeat_data(:, end));
        
        %Normalize by initial count (% remaining)
        percent_remaining(r_idx, t_idx) = 100 * (final_mean / initial_count);
    end
end

%Plotting heatmap
figure;
h = heatmap(totals, ratios, percent_remaining);
h.Title = 'Nucleus Pulposus % Remaining';
h.XLabel = 'Total Macrophage Count';
h.YLabel = 'M2:M1 Ratio';
colormap parula;  % optional color scheme
colorbar;

%% THIS IS MY FIGURE 2
%Fits story to evaluate independent effects of each parameter first
%This will yield 4 graphs (2 each for NP over time, and 2 each for endpoint % remaining)
%One of each will be varying ratio while holding total constant and vice versa

% time conversion to simulation days
dt = 60;
time_days = (0:1000) * dt / 1440;

ratios = [0, 0.2, 0.4, 0.5, 0.6, 0.8, 1];
totals = [100, 200, 300, 400];
num_repeats = 3;
initial_count = 4009;

base_purple = [59 53 97] / 255; %this is my purple color I like

%GRAPH 1 - NP over time (vary ratio, total = 300) 

figure;
hold on
t_idx = 3; % 300 macrophages

for r_idx = 1:length(ratios)
    
    start_row = (t_idx-1)*length(ratios)*num_repeats + (r_idx-1)*num_repeats + 1;
    end_row = start_row + num_repeats - 1;
    
    data = secondfile(start_row:end_row, :);
    mean_curve = mean(data,1);
    
    % gradient color (lighter for lower ratios)
    alpha = r_idx / length(ratios);
    color = base_purple * alpha + (1-alpha)*[1 1 1];
    
    plot(time_days, mean_curve, 'Color', color, 'LineWidth', 2);
end

xlabel('Simulation Time (days)')
ylabel('Nucleus Pulposus Population')
title('NP Population vs Time (Total Macrophages = 300)')
legend(string(ratios)+ " M2:M1", 'Location','northwest')
xlim([0 max(time_days)])

%GRAPH 2 - Endpoint NP % (vary ratio, total = 300) 

figure;
hold on

t_idx = 3; % total = 300

means = zeros(1,length(ratios));
stds = zeros(1,length(ratios));

all_points = cell(1,length(ratios));

for r_idx = 1:length(ratios)
    
    start_row = (t_idx-1)*length(ratios)*num_repeats + (r_idx-1)*num_repeats + 1;
    end_row = start_row + num_repeats - 1;
    
    data = secondfile(start_row:end_row, :);
    final_vals = 100 * data(:,end) / initial_count;
    
    means(r_idx) = mean(final_vals);
    stds(r_idx) = std(final_vals);
    
    all_points{r_idx} = final_vals;
end

%bar plot
b = bar(means);
b.FaceColor = 'flat';

for i = 1:length(ratios)
    alpha = i / length(ratios);
    b.CData(i,:) = base_purple * alpha + (1-alpha)*[1 1 1];
end

%error bars
errorbar(1:length(ratios), means, stds, 'k.', 'LineWidth',1.5)

xticks(1:length(ratios))
xticklabels(string(ratios))
xlabel('M2:M1 Ratio')
ylabel('NP % Remaining')
title('Endpoint NP % vs M2:M1 Ratio (Total = 300)')

%statistics
values = [];
group = [];

for r_idx = 1:length(ratios)
    values = [values; all_points{r_idx}];
    group = [group; r_idx*ones(num_repeats,1)];
end

[p,~,stats] = anova1(values, group, 'off');
results = multcompare(stats, 'Display','off');

base_gap = max(stds) * 0.15;
lane_offset = max(stds) * 0.5;

lane_counter = containers.Map('KeyType','char','ValueType','double');

for i = 1:size(results,1)

    g1 = results(i,1);
    g2 = results(i,2);
    p_val = results(i,6);

    is_adjacent = abs(g1 - g2) == 1;
    is_vs_zero  = (g1 == 1 || g2 == 1);

    if ~(is_adjacent || is_vs_zero)
        continue
    end

    if p_val < 0.001
        stars = '***';
    elseif p_val < 0.01
        stars = '**';
    elseif p_val < 0.05
        stars = '*';
    else
        continue
    end

    %this helps with positioning of the significance labeling
    anchor = max(means([g1 g2]) + stds([g1 g2]));

    is_full_span = (g1 == 1 && g2 == length(ratios)) || ...
                   (g2 == 1 && g1 == length(ratios));

    if is_full_span
        anchor = means(end) + stds(end);
    end
    
    %this tiers the lines
    if is_vs_zero && ~is_adjacent
        tier = 2;
    else
        tier = 1;
    end

    tier_key = num2str(tier);

    %creating lanes to separate tiers
    if isKey(lane_counter, tier_key)
        lane_counter(tier_key) = lane_counter(tier_key) + 1;
    else
        lane_counter(tier_key) = 0;
    end

    lane = lane_counter(tier_key);

    %final y position of stats on graph
    if is_full_span
        % keep 0 vs 1 tight to last bar
        y = anchor + base_gap + (lane * lane_offset * 0.3);
    else
       y = anchor + base_gap + (lane * lane_offset * 0.3) + max(stds)*0.007;
    end

    if (g1 == 1 && g2 == 3) || (g1 == 3 && g2 == 1)
        y = y + max(stds) * 0.6;
    end

    if (g1 == 1 && g2 == 4) || (g1 == 4 && g2 == 1)
        y = y - max(stds) * 0.6;
    end

    %actually drawing the lines and stars
    plot([g1 g2], [y y], 'k', 'LineWidth',1.2)

    text(mean([g1 g2]), y + 0.15, stars, ...
        'HorizontalAlignment','center', ...
        'FontSize',18)
end

%starting the bar at 0
ylim([0, max(means+stds) + max(stds)*2])



%GRAPH 3 - NP over time (vary total, ratio = 0.5)

figure;
hold on

r_idx = 4; % ratio = 0.5

for t_idx = 1:length(totals)
    
    start_row = (t_idx-1)*length(ratios)*num_repeats + (r_idx-1)*num_repeats + 1;
    end_row = start_row + num_repeats - 1;
    
    data = secondfile(start_row:end_row, :);
    mean_curve = mean(data,1);
    
    alpha = t_idx / length(totals);
    color = base_purple * alpha + (1-alpha)*[1 1 1];
    
    plot(time_days, mean_curve, 'Color', color, 'LineWidth', 2);
end

xlabel('Simulation Time (days)')
ylabel('Nucleus Pulposus Population')
title('NP Population vs Time (M2:M1 = 0.5)')
legend(string(totals)+ " Total Macrophages", 'Location','northeast')
xlim([0 max(time_days)])

%GRAPH 4 - Endpoint NP % (vary total, ratio = 0.5)

figure;
hold on

r_idx = 4; % ratio = 0.5

%initializing vectors
means = zeros(1,length(totals));
stds  = zeros(1,length(totals));

values = [];
group  = [];

for t_idx = 1:length(totals)

    start_row = (t_idx-1)*length(ratios)*num_repeats + (r_idx-1)*num_repeats + 1;
    end_row   = start_row + num_repeats - 1;

    data = secondfile(start_row:end_row, :);
    final_vals = 100 * data(:,end) / initial_count;

    means(t_idx) = mean(final_vals);
    stds(t_idx)  = std(final_vals);

    values = [values; final_vals];
    group  = [group; t_idx*ones(num_repeats,1)];
end

%creating bar graph
b = bar(means, 'FaceColor','flat');

base_purple = [0.25 0.2 0.45]; 

for i = 1:length(totals)
    alpha = i / length(totals);
    b.CData(i,:) = base_purple * alpha + (1-alpha)*[1 1 1];
end

errorbar(1:length(totals), means, stds, 'k.', 'LineWidth', 1.5)

xticks(1:length(totals))
xticklabels(string(totals))
xlabel('Total Macrophage Count')
ylabel('NP % Remaining')
title('Endpoint NP % vs Total Macrophages (M2:M1 = 0.5)')

%doing the stastistics
[p,~,stats] = anova1(values, group, 'off');
results = multcompare(stats, 'Display','off');


%placing significance stars
y_base = max(means) + max(stds)*1.5;
y_step = max(stds)*0.8;

sig_map = containers.Map('KeyType','char','ValueType','double');

for i = 1:size(results,1)

    g1 = results(i,1);
    g2 = results(i,2);
    p_val = results(i,6);

    if p_val >= 0.05
        continue
    elseif p_val < 0.001
        stars = '***';
    elseif p_val < 0.01
        stars = '**';
    else
        stars = '*';
    end

    %position just above taller bar
    y = max(means([g1 g2]) + stds([g1 g2])) + 1;
    base_gap = max(stds) * 0.6;   %slight vertical lift above bar

    %anchor to taller bar in comparison
    anchor = max(means(g1), means(g2));

    dist = abs(g1 - g2);

        %structured spacing tiers
        if dist == 3
            tier = 3;   %100 vs 400 (highest)
        elseif dist == 2
            tier = 2;   %100 vs 300, 200 vs 400
        elseif dist == 1
            tier = 1;   %adjacent comparisons
        else
            tier = 1;
        end

        %lane system to prevent overlap within same tier
        if ~exist('lane_counter','var')
            lane_counter = containers.Map('KeyType','char','ValueType','double');
        end
        
        key = sprintf('%d-%d', g1, g2);
        
        if isKey(lane_counter, num2str(tier))
            lane_counter(num2str(tier)) = lane_counter(num2str(tier)) + 1;
        else
            lane_counter(num2str(tier)) = 1;
        end
        
        lane = lane_counter(num2str(tier));
        
        %final y position of lines and stars
        y = anchor + base_gap + (tier * 1.6) + (lane * 0.8 ...
            );
    plot([g1 g2], [y y], 'k', 'LineWidth',1.2)
    
    %plotting stars
    text(mean([g1 g2]), y + 0.3, stars, ...
        'HorizontalAlignment','center', ...
        'FontSize',18)
end

% ylim([min(means)-5, y_base + y_step*(length(sig_map)+2)])
ylim([0, max(means+stds) + max(stds)*2])
box off
hold off


%% %% Sensitivity Analysis for % Remaining Nucleus Pulposus
ratios = [0, 0.2, 0.4, 0.5, 0.6, 0.8, 1];
totals = [100, 200, 300, 400];
num_repeats = 3;
num_ratios = length(ratios);
num_totals = length(totals);
initial_count = 4009;

percent_remaining = nan(num_ratios, num_totals);
percent_remaining_trials = nan(num_ratios, num_totals, num_repeats);

for t_idx = 1:num_totals
    for r_idx = 1:num_ratios
         % Row indices in secondfile for this total x ratio
        start_row = (t_idx-1)*num_ratios*num_repeats + (r_idx-1)*num_repeats + 1;
        end_row = start_row + num_repeats - 1;

        repeat_data = secondfile(start_row:end_row, :);
        final_vals = repeat_data(:, end);
        percent_remaining_trials(r_idx, t_idx, :) = 100 * (final_vals / initial_count);
        percent_remaining(r_idx, t_idx) = mean(percent_remaining_trials(r_idx, t_idx, :));
    end
end


%getting in vectors to run the ANOVA
ratio_vec = [];
total_vec = [];
NP_vec = [];

for t_idx = 1:num_totals
    for r_idx = 1:num_ratios
        for rep = 1:num_repeats
            ratio_vec(end+1) = ratios(r_idx);
            total_vec(end+1) = totals(t_idx);
            NP_vec(end+1) = percent_remaining_trials(r_idx, t_idx, rep);
        end
    end
end


%ANOVA
[p, tbl, stats] = anovan(NP_vec, {ratio_vec, total_vec}, ...
    'model','interaction', ...
    'varnames',{'M2_M1_ratio','Macrophage_count'});
%matlab makes its own table
disp('ANOVA Table:');
disp(tbl);

%% FIGURE 3B
base_purple = [59 53 97] / 255;

figure; hold on

for i = 1:num_totals
    
    % create gradient color (light → dark)
    alpha = i / num_totals;
    color_i = base_purple * alpha + (1-alpha)*[1 1 1];
    
    errorbar(ratios, meanNP(:,i), stdNP(:,i), ...
        '-o', ...
        'Color', color_i, ...
        'LineWidth', 2, ...
        'MarkerSize', 6, ...
        'CapSize', 8, ...
        'DisplayName', sprintf('%d Total Macrophages', totals(i)));
end
set(gca, 'FontSize', 14)
xlabel('M2:M1 Ratio')
ylabel('Endpoint % NP Remaining')
legend('Location','southeast')
lgd.ItemTokenSize = [25, 20];




