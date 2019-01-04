function crunch_rouder2004data(varargin)

% Inputs
input_method = varargin{1};
switch input_method
    case 1
        filename = varargin{2};
        load(filename);
    case 2
        Cd = varargin{2};
        %         Ch = varargin{3};
        %         Cs = varargin{4};
end

%% Figure 3a
count = 0;
frac_correct = cell(3,1);
for nstim = [13 20 30]
    count = count + 1;
    frac_correct{count} = zeros(nstim, 1);
    experiment_index_set = Cd{8}==nstim;
    for stim_size = 1:nstim
        index_set = experiment_index_set & Cd{4}==stim_size-1;
        nr_all = sum(index_set);
        nr_correct = sum(Cd{5}(index_set)==stim_size-1);
        frac_correct{count}(stim_size) = nr_correct / nr_all;        
    end
end
figure,
plot(1:13, frac_correct{1}, 'ko-', 'MarkerSize', 10);
hold on;
plot(1:20, frac_correct{2}, 'ks-', 'MarkerSize', 10, ...
    'MarkerFaceColor', [0.7 0.7 0.7]);
plot(1:30, frac_correct{3}, 'kd-', 'MarkerSize', 10, ...
    'MarkerFaceColor', [0 0 0]);
xlim([0 31]);
ylim([0 1]);
xlabel('Stimulus (length)');
ylabel('Fraction of correct at first guess');
set(gca, 'FontSize', 20);
legend('N=13', 'N=20', 'N=30', 'Location', 'Best');

%% Figure 4 - Efficiency figures
cost_of_error_type = 3;
subjects = {'Monique', 'Richard'};
nr_subjets = length(subjects);
smin = -1e2;
smax = -1e-8;
for nstim = [13 20 30]
    D_emp = zeros(nr_subjets, 1);
    empirical_info_rate = zeros(nr_subjets, 1);
    efficiency = zeros(nr_subjets, 1);
    for k = 1:2
        name = subjects{k};
        DR_simulated = findRate([smin smax 1000], ...
            ['N=' num2str(nstim) ' absolute identification'], ...
            'cost_of_error_type', cost_of_error_type);
        xlim([0 max(1, max(DR_simulated(:,1)))]);
        ylim([0 6]);
        hold on;
        [D_emp(k), empirical_info_rate(k), efficiency(k)] = ...
            draw_efficiency(nstim, name, DR_simulated, Cd);
    end
    % Both efficiencies plotted    
    findRate([smin smax 1000], ...
        ['N=' num2str(nstim) ' absolute identification'], ...
        'cost_of_error_type', cost_of_error_type);
    xlim([0 max(1, max(DR_simulated(:,1)))]);
    ylim([0 6]);
    hold on;
    for k = 1:2
        plot(D_emp(k), empirical_info_rate(k), 'k.', 'MarkerSize', 20);
        text(D_emp(k) + 0.025, empirical_info_rate(k) + 0.25, ...
            [num2str(efficiency(k), '%1.3f')], ...
            'FontSize', 15);
    end

end

% Figure 4* - everyone's efficieny
subjects = {'everyone'};
for nstim = [13 20 30]
    name = subjects{1};
    DR_simulated = findRate([smin smax 1000], ...
        ['N=' num2str(nstim) ' absolute identification'], ...
        'cost_of_error_type', cost_of_error_type);
    xlim([0 max(1, max(DR_simulated(:,1)))]);
    ylim([0 6]);
    hold on;
    draw_efficiency(nstim, name, DR_simulated, Cd);
end

%% Figure 5a
human_capacity_limits = [3.01, 3.25, 3.19];
count = 0;
frac_correct = cell(3,1);
optimal_Qyx_frac_correct = cell(3,1);
for nstim = [13 20 30]
    % empirical data
    count = count + 1;
    frac_correct{count} = zeros(nstim, 1);
    experiment_index_set = Cd{8}==nstim;
    for stim_size = 1:nstim
        index_set = experiment_index_set & Cd{4}==stim_size-1;
        nr_all = sum(index_set);
        nr_correct = sum(Cd{5}(index_set)==stim_size-1);
        frac_correct{count}(stim_size) = nr_correct / nr_all;        
    end
    
    % optimal Qyx
    cap_lim = human_capacity_limits(count);
    [~, Qyx_optimal] = findRate([-1e1 -1e-4 1000], ...
        ['N=' num2str(nstim) ' absolute identification'], ...
        'capacity_limit', cap_lim);
    optimal_Qyx_frac_correct{count} = Qyx_optimal(1:nstim+1:end);
    
end
figure,
% empirical data
plot(1:13, frac_correct{1}, 'ko', 'MarkerSize', 10);
hold on;
plot(1:20, frac_correct{2}, 'ks', 'MarkerSize', 10, ...
    'MarkerFaceColor', [0.7 0.7 0.7]);
plot(1:30, frac_correct{3}, 'kd', 'MarkerSize', 10, ...
    'MarkerFaceColor', [0 0 0]);
xlim([0 31]);
ylim([0 1]);
xlabel('Stimulus (length)');
ylabel('Fraction of correct at first guess');
set(gca, 'FontSize', 20);
% optimal Qyx
plot(1:13, optimal_Qyx_frac_correct{1}, 'r');
plot(1:20, optimal_Qyx_frac_correct{2}, 'r');
plot(1:30, optimal_Qyx_frac_correct{3}, 'r');
legend('N=13', 'N=20', 'N=30', 'Model L_1', 'Location', 'Best');

%% Figure 5b
human_capacity_limits = [3.01, 3.25, 3.19];
count = 0;
frac_correct = cell(3,1);
optimal_Qyx_frac_correct = cell(3,1);
for nstim = [13 20 30]
    % empirical data
    count = count + 1;
    frac_correct{count} = zeros(nstim, 1);
    experiment_index_set = Cd{8}==nstim;
    for stim_size = 1:nstim
        index_set = experiment_index_set & Cd{4}==stim_size-1;
        nr_all = sum(index_set);
        nr_correct = sum(Cd{5}(index_set)==stim_size-1);
        frac_correct{count}(stim_size) = nr_correct / nr_all;        
    end
    
    % optimal Qyx
    cap_lim = human_capacity_limits(count);
    [~, Qyx_optimal] = findRate([-1e1 -1e-4 1000], ...
        ['N=' num2str(nstim) ' absolute identification'], ...
        'capacity_limit', cap_lim, 'cost_of_error_type', 2);
    optimal_Qyx_frac_correct{count} = Qyx_optimal(1:nstim+1:end);
    
end
figure,
% empirical data
plot(1:13, frac_correct{1}, 'ko', 'MarkerSize', 10);
hold on;
plot(1:20, frac_correct{2}, 'ks', 'MarkerSize', 10, ...
    'MarkerFaceColor', [0.7 0.7 0.7]);
plot(1:30, frac_correct{3}, 'kd', 'MarkerSize', 10, ...
    'MarkerFaceColor', [0 0 0]);
xlim([0 31]);
ylim([0 1]);
xlabel('Stimulus (length)');
ylabel('Fraction of correct at first guess');
set(gca, 'FontSize', 20);
legend('N=13', 'N=20', 'N=30', 'Location', 'Best');
% optimal Qyx
plot(1:13, optimal_Qyx_frac_correct{1}, 'r');
plot(1:20, optimal_Qyx_frac_correct{2}, 'r');
plot(1:30, optimal_Qyx_frac_correct{3}, 'r');
legend('N=13', 'N=20', 'N=30', 'Model L_2', 'Location', 'Best');


%% Figure 5c
cost_of_error_type = 1;
human_capacity_limits = [3.01, 3.25, 3.19];
count = 0;
frac_correct = cell(3,1);
optimal_Qyx_frac_correct = cell(3,1);
for nstim = [13 20 30]
    % empirical data
    count = count + 1;
    frac_correct{count} = zeros(nstim, 1);
    experiment_index_set = Cd{8}==nstim;
    for stim_size = 1:nstim
        index_set = experiment_index_set & Cd{4}==stim_size-1;
        nr_all = sum(index_set);
        nr_correct = sum(Cd{5}(index_set)==stim_size-1);
        frac_correct{count}(stim_size) = nr_correct / nr_all;        
    end
    
    % optimal Qyx
    cap_lim = human_capacity_limits(count);
    [~, Qyx_optimal] = findRate([-1e2 -1e-3 1000], ...
        ['N=' num2str(nstim) ' absolute identification'], ...
        'capacity_limit', cap_lim, ...
        'cost_of_error_type', cost_of_error_type);
    optimal_Qyx_frac_correct{count} = Qyx_optimal(1:nstim+1:end);
    
end
figure,
% empirical data
plot(1:13, frac_correct{1}, 'ko', 'MarkerSize', 10);
hold on;
plot(1:20, frac_correct{2}, 'ks', 'MarkerSize', 10, ...
    'MarkerFaceColor', [0.7 0.7 0.7]);
plot(1:30, frac_correct{3}, 'kd', 'MarkerSize', 10, ...
    'MarkerFaceColor', [0 0 0]);
xlim([0 31]);
ylim([0 1]);
xlabel('Stimulus (length)');
ylabel('Fraction of correct at first guess');
set(gca, 'FontSize', 20);
legend('N=13', 'N=20', 'N=30', 'Location', 'Best');
% optimal Qyx
plot(1:13, optimal_Qyx_frac_correct{1}, 'r');
plot(1:20, optimal_Qyx_frac_correct{2}, 'r');
plot(1:30, optimal_Qyx_frac_correct{3}, 'r');
legend('N=13', 'N=20', 'N=30', 'Model L_3', 'Location', 'Best');


%% Human capacity limits
for nstim = [13 20 30]
    DR_simulated = findRate([-1e1 -1e-4 1000], ...
        ['N=' num2str(nstim) ' absolute identification']);
    [~, empirical_info_rate] = draw_efficiency(nstim, 'everyone', ...
        DR_simulated, Cd);
    empirical_info_rate
end

%% Figure 7
f = @(x,n) 0.5*( tan((2*(x-1)/(n+1)-1)*pi/4) + 1 );
reza = @(x,n) f(x,n)/f(n,n);
n = 30;
cost_matrix = zeros(n);
scale = 4;
for i=1:n*scale
    for j=1:n*scale
        cost_matrix(i,j) = abs(reza(j/scale,n)-reza(i/scale,n));
    end
end
figure,
contourf(cost_matrix, 40);
colormap(jet(256));
colorbar('eastoutside');
set(gca,'XTick',scale:5*scale:n*scale,'XTickLabel',1:5:n)
set(gca,'YTick',scale:5*scale:n*scale,'YTickLabel',1:5:n)
set(gca, 'FontSize', 20);
xlabel('Stimulus');
ylabel('Response');
title('$L(x,y) = |r(y)-r(x)|$', 'Interpreter', 'Latex');

%% Figure 6b
figure,
plot((1:13)/13, reza(1:13,13), 'ko-');
hold on;
plot((1:20)/20, reza(1:20,20), 'ks-');
plot((1:30)/30, reza(1:30,30), 'kd-');
set(gca, 'FontSize', 20);
xlabel('Fraction of stimulus');
ylabel('Anchor')
legend('N=13', 'N=20', 'N=30', 'Location', 'Best');

end

