function [DR, Qyx_out, capacity_limit]= findRate(varargin)

% Inputs
smin = varargin{1}(1);
smax = varargin{1}(2);
nr_s = varargin{1}(3);
experiment_type = varargin{2};
capacity_limit = inf;
cost_of_error_type = 1;
% ----------------------------------------------------------------------
i=3;
while i<=length(varargin),
    switch varargin{i},
        case 'capacity_limit',
            capacity_limit = varargin{i+1};
        case 'cost_of_error_type'
            cost_of_error_type = varargin{i+1};
        otherwise,
            display(varargin{i});
            error('Unexpected inputs!!!');
    end
    i = i+2;
end
% ----------------------------------------------------------------------
% Example for dice
switch experiment_type
    case 'N=30 absolute identification'
        n = 30;                         % # of sides        
        p = ones(n,1) * 1/n;            % fair dice
    case 'N=20 absolute identification'
        n = 20;                         % # of sides        
        p = ones(n,1) * 1/n;            % fair dice
    case 'N=13 absolute identification'
        n = 13;                         % # of sides        
        p = ones(n,1) * 1/n;            % fair dice
    case 'N=12 color identification'
        n = 12;                         % # of sides        
        p = ones(n,1) * 1/n;            % fair dice        
    case 'dice'
        n = 6;                          % # of sides
        p = ones(n,1) * 1/n;            % fair dice
end

switch cost_of_error_type
    case 1
        cost_matrix = ones(n);          % L(x, y) = 1 if x ~= y
        cost_matrix(1:n+1:end) = 0;     % L(x, y) = 0 if x  = y
    case 2
        % L(x, y) = |y-x|
        cost_matrix = zeros(n);
        for i=1:n
            for j=1:n
                cost_matrix(i,j) = abs(j-i); 
            end
        end
    case 3
        % L(x, y) = |reza(y)-reza(x)|
        f = @(x,n) 0.5*( tan((2*(x-1)/(n+1)-1)*pi/4) + 1 );
        reza = @(x,n) f(x,n)/f(n,n);
        cost_matrix = zeros(n);
        for i=1:n
            for j=1:n
                cost_matrix(i,j) = abs(reza(j,n)-reza(i,n)); 
            end
        end
end

Qyx_init = eye(n);              % initial guess

s_vec = (-1) * logspace(log10(-smax), log10(-smin), nr_s);
DR = zeros(nr_s, 2);
Qyx = cell(nr_s);
for k = 1:nr_s
    current_s_val = s_vec(k);
    [DR(k, 1), DR(k, 2), Qyx{k}] = calculate_rate_and_distortion(...
        current_s_val, cost_matrix, Qyx_init, p);
end    

% Plot rate-distortion curve
figure,
plot(DR(:, 1), DR(:, 2), 'o');
xlabel('Distortion');
ylabel('Rate');
set(gca, 'FontSize',20);

switch experiment_type
    case 'dice'
        % Plot theoretical curve
        hold on;
        Dmin = min(DR(:,1));
        Dmax = max(DR(:,1));
        nr_D = 100;
        D_vec = linspace(Dmin, Dmax, nr_D);
        R_theory = log2(6-6*D_vec) + D_vec .* log2( D_vec ./ (5-5*D_vec) );
        plot(D_vec, R_theory, 'r');        
        % Add legend
        legend('Blahut-Arimoto algorithm', ...
            'Theoretical', 'Location', 'Best');

    otherwise
        % Add legend
        legend('Blahut-Arimoto algorithm', 'Location', 'Best');
        fprintf('FO: No theoretical rate-distortion curve.\n');

end

% Add title
title(experiment_type, 'FontSize', 20);

% Find optimal Qyx if a capacity is given
if capacity_limit < inf
    [~, idx] = min(abs(DR(:,2)-capacity_limit));
    Qyx_out = Qyx{idx};
end

end

%%
function [D, R, Qyx] = calculate_rate_and_distortion(s, cost_matrix, Qyx, p)
display_result = 0;
nr_iter = 100;

A = exp(s*cost_matrix);         % this term will frequently be used

sz = size(Qyx);
nr_rows = sz(1);
nr_columns = sz(2);
q = zeros(nr_columns,1);
for iter_count = 1:nr_iter
    % Step 1: Calculate q(y_j) from Q(y_j|x_i)
    %         Qyx(i,j) = Q(y_j|x_i)    
    for j = 1:nr_columns
        q(j) = sum(Qyx(:,j) .* p);
    end
    
    % Step 2: Calculate Q(y_j|x_i) from q(y_j)
    for j = 1:nr_columns
        Qyx(:,j) = q(j) * A(:,j);
    end
    rowsum = sum(Qyx, 2);
    for i = 1:nr_rows
        Qyx(i,:) = Qyx(i,:) / rowsum(i);
    end
    
    % Print-out iteration result
    if display_result
        fprintf('Iter %3i: q = [', iter_count);
        fprintf('%1.6f ', q);
        fprintf('\n');
    end
end

% Calculate lambda (Lagrance multiplier)
lambda = zeros(nr_rows,1);
for i = 1:nr_rows
    lambda(i) = 1 / sum( q' .* A(i,:) );
end

% Calculate D-distortion and R-rate, as values dependent on s
D = 0;
for i = 1:nr_rows
    for j = 1:nr_columns
        D = D + lambda(i) * p(i) * q(j) * A(i,j) * cost_matrix(i,j);
    end
end
R = s*D;
for i = 1:nr_rows
    R = R + sum( p(i) .* log(lambda(i)) );
end
R = R * log2(exp(1));   % this script computes the rate in nats

end

