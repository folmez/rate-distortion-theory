function crunch_nosofsky1987data

load('nosofsky1987data.mat');
nstim = length(nosofsky1987_Cd);

%% Plot fraction of correct
frac_correct = nosofsky1987_Cd(1:nstim+1:end)' ./ sum(nosofsky1987_Cd, 2);
figure,
plot(1:nstim, frac_correct, 'ko-')
ylim([0 1]);
set(gca, 'FontSize', 20);
xlabel('color');
ylabel('fraction of correct guess');
title(['N=' num2str(nstim) ' color identification'], 'FontSize', 20);

%% Rate-distortion curve and efficiency
DR_simulated = findRate([-1e1 -1e-4 1000], ...
    ['N=' num2str(nstim) ' color identification']);
xlim([0 1]);
ylim([0 4]);
hold on;

% max distortion
D_max = max(DR_simulated(:,1));

% calculate D_emp, empirical distortion
nr_all = sum(sum(nosofsky1987_Cd));
nr_fail = nr_all - sum(nosofsky1987_Cd(1:nstim+1:end));
D_emp = nr_fail / nr_all;
plot([D_emp D_emp], [0 6], 'k');

% D_star
% empirical P(y_j|x_i) = pyx(i,j)
pyx_emp = zeros(nstim);
for i = 1:nstim
    pyx_emp(i, :) = nosofsky1987_Cd(i, :) / ...
        sum(nosofsky1987_Cd(i, :));
end
% empirical i/o probability
px_emp = sum(nosofsky1987_Cd, 2);
py_emp = sum(nosofsky1987_Cd, 1)';
px_emp = px_emp / sum(px_emp);
py_emp = py_emp / sum(py_emp);
% calculate empirical info rate according to Eq. (2)
empirical_info_rate = 0;
for i = 1:nstim
    for j = 1:nstim
        if pyx_emp(i,j)~=0
            empirical_info_rate = empirical_info_rate + ...
                pyx_emp(i, j) * px_emp(i) * log2(pyx_emp(i,j)/py_emp(j));
        end
    end
end
[~, idx] = min(abs(DR_simulated(:,2)-empirical_info_rate));
R_star = DR_simulated(idx, 2);
plot([0 1], [R_star R_star], 'k');
D_star = DR_simulated(idx, 1);
plot([D_star D_star], [0 6], 'k--');
efficiency = (D_emp - D_max) / (D_star - D_max);
plot(D_emp, empirical_info_rate, 'k.', 'MarkerSize', 20);
text(D_emp + 0.025, empirical_info_rate + 0.25, ...
    ['efficiency = ' num2str(efficiency, '%1.3f')], ...
    'FontSize', 15);
    

%% 
human_capacity_limit = empirical_info_rate;
frac_correct = zeros(nstim, 1);
for stim_size = 1:nstim
    nr_all = sum(nosofsky1987_Cd(stim_size, :));
    nr_correct = nosofsky1987_Cd(stim_size, stim_size);
    frac_correct(stim_size) = nr_correct / nr_all;
end

% optimal Qyx
[~, Qyx_optimal] = findRate([-1e1 -1e-4 1000], ...
    ['N=' num2str(nstim) ' color identification'], ...
    'capacity_limit', human_capacity_limit);
optimal_Qyx_frac_correct = Qyx_optimal(1:nstim+1:end);

figure,
% empirical data
plot(1:nstim, frac_correct, 'ko', 'MarkerSize', 10);
hold on;
xlim([0 nstim+1]);
ylim([0 1]);
xlabel('Stimulus (color)');
ylabel('Fraction of correct at first guess');
set(gca, 'FontSize', 20);
% optimal Qyx
plot(1:nstim, optimal_Qyx_frac_correct, 'r');
legend('N=12', 'Model L_1', 'Location', 'Best');


end