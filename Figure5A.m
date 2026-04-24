%The purpose of this file is to create Figure 5A
%This graph plots the means of each treatment group as well as their corresponding bootstrapped 95% confidence intervals
%The untreated line is plotted for the total 60,000 minutes. the treatments are initiated at 15,000 minutes in
clear; close all;

%defining time
dt = 60;
time_minutes = (0:1000) * dt;
time_days = time_minutes / 1440;

start_idx = 251;

t_post = time_days(start_idx:end);

%defining my colors
color_untreated = [59 53 97] / 255;
color_sealant = [219 43 57] / 255;
color_cytokine = [121 181 153] / 255;
color_celltherapy = [0 141 213] / 255;

%untreated control data
load('kyanna30trialsdata.mat');

untreatedcontrol_NP = nSTORAGE(:,:);
untreatedcontrol_NP_data = untreatedcontrol_NP(4:33,:);

untreatedcontrol_NP_mean_curve = mean(untreatedcontrol_NP_data,1);
B = 1000;
n_time = size(untreatedcontrol_NP_data,2);

untreatedcontrol_ci_lower = zeros(1,n_time);
untreatedcontrol_ci_upper = zeros(1,n_time);

for t = 1:n_time
    ci = bootci(B,@mean,untreatedcontrol_NP_data(:,t));
    untreatedcontrol_ci_lower(t) = ci(1);
    untreatedcontrol_ci_upper(t) = ci(2);
end


%HIGH DENSITY SEALANT
load('kyannahighdensitywound2.mat');

sealant_NP_portion = nSTORAGE(1:30,:);
sealant_NP_combined = sealant_NP_portion;

sealant_NP_mean_curve = mean(sealant_NP_combined,1);

n_time = size(sealant_NP_combined,2);

sealant_ci_lower = zeros(1,n_time);
sealant_ci_upper = zeros(1,n_time);

for t = 1:n_time
    ci = bootci(B,@mean,sealant_NP_combined(:,t));
    sealant_ci_lower(t) = ci(1);
    sealant_ci_upper(t) = ci(2);
end

%M2 CELL THERAPY
load('kyannam2celltherapy.mat');

celltherapy_NP_portion = nSTORAGE(1:30,:);
celltherapy_NP_combined = celltherapy_NP_portion;

celltherapy_NP_mean_curve = mean(celltherapy_NP_combined,1);

n_time = size(celltherapy_NP_combined,2);

celltherapy_ci_lower = zeros(1,n_time);
celltherapy_ci_upper = zeros(1,n_time);

for t = 1:n_time
    ci = bootci(B,@mean,celltherapy_NP_combined(:,t));
    celltherapy_ci_lower(t) = ci(1);
    celltherapy_ci_upper(t) = ci(2);
end

%CYTOKINE LOADED GEL
load('kyannasoftgel.mat');

cytokine_NP_portion = nSTORAGE(1:30,:);
cytokine_NP_combined = cytokine_NP_portion;

cytokine_NP_mean_curve = mean(cytokine_NP_combined,1);

n_time = size(cytokine_NP_combined,2);

cytokine_ci_lower = zeros(1,n_time);
cytokine_ci_upper = zeros(1,n_time);

for t = 1:n_time
    ci = bootci(B,@mean,cytokine_NP_combined(:,t));
    cytokine_ci_lower(t) = ci(1);
    cytokine_ci_upper(t) = ci(2);
end

%PLOTTING
figure;
set(gcf,'Position',[100 100 950 600])
hold on

%untreated
fill([time_days fliplr(time_days)], ...
     [untreatedcontrol_ci_upper fliplr(untreatedcontrol_ci_lower)], ...
     color_untreated,'FaceAlpha',0.15,'EdgeColor','none');

h1 = plot(time_days, untreatedcontrol_NP_mean_curve, ...
     'Color',color_untreated,'LineWidth',2);


%high density sealant
fill([t_post fliplr(t_post)], ...
     [sealant_ci_upper fliplr(sealant_ci_lower)], ...
     color_sealant,'FaceAlpha',0.15,'EdgeColor','none');

h2 = plot(t_post, sealant_NP_mean_curve, ...
     'Color',color_sealant,'LineWidth',2);

%M2 cell therapy
fill([t_post fliplr(t_post)], ...
     [celltherapy_ci_upper fliplr(celltherapy_ci_lower)], ...
     color_celltherapy,'FaceAlpha',0.15,'EdgeColor','none');

h3 = plot(t_post, celltherapy_NP_mean_curve, ...
     'Color',color_celltherapy,'LineWidth',2);

%cytokine loaded gel
fill([t_post fliplr(t_post)], ...
     [cytokine_ci_upper fliplr(cytokine_ci_lower)], ...
     color_cytokine,'FaceAlpha',0.15,'EdgeColor','none');

h4 = plot(t_post, cytokine_NP_mean_curve, ...
     'Color',color_cytokine,'LineWidth',2);

%Formatting the plot

xline(time_days(start_idx),'--k','Intervention', 'FontSize', 14);

xlabel('Simulation Time (days)','FontSize', 14)
ylabel('Mean Nucleus Pulposus Population', 'FontSize', 14)
%title('Mean Nucleus Pulposus Population Over Simulated Time')

xlim([0 max(time_days)])

lgd = legend([h1 h2 h3 h4], ...
    {'Untreated','High Density Sealant', ...
     'M2 Cell Therapy','Cytokine Loaded Gel'}, ...
    'Location','northeast')
lgd.FontSize = 12;
