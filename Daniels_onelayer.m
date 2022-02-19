
clear; clc; close all;


%% Set parameters
rng(117)
num_bar = 15;
force_mean = (num_bar)*200; % External force kN
Statistical_data = ones(num_bar); % Area information, mm2
Demand_cov = 0.35;
Demand_mean = 400; % Newton
MCS_N = 1e6;

%% Estimate reliability index

% Estimate failure probability for a single bar
% Estimate Beta
tmp_sig = Demand_mean*Demand_cov;
tmp_beta = (Demand_mean - force_mean/num_bar)/tmp_sig;
%(400-3000/15)/(400*0.35)
% Estimate failure probability    
fail_prob = normcdf(-tmp_beta); % Failure probability

total_reli_prob = [];
for ii = 1: num_bar
    
    tmp = fail_prob^(ii)*(1-fail_prob)^(num_bar-ii);
    total_reli_prob = [total_reli_prob; tmp];
    
end

Relia_index = [];
for ii =1:length(total_reli_prob)
    tmp = total_reli_prob(ii);
    Relia_index = [Relia_index; -norminv(tmp)];
end

%% Estimate redundancy index
total_redu_prob = [];

for ii = 1: num_bar
    
    % when all elements fail
    if ii == num_bar

        % The failure cases where all elements are not existed
        non_fail = 0;
        fail = MCS_N;
    else

        % Failure cases where some elements are survived
        % Load for each survived element
        load_init = ones([1,num_bar-ii])*force_mean/(num_bar-ii);
        tmp_total_C = 1:1:length(Statistical_data)-ii;
        
        % MCS is carried out
        non_fail = 0; fail = 0; %tmp_prob_history = zeros(MCS_N,1);
        for kk = 1 : MCS_N

            % Randomly generate capacities
            capacity = normrnd(Demand_mean, Demand_mean*Demand_cov, [1,length(load_init)]);

            % Check redundancy
            redundancy_true = true;       % While loop indicator
            demand = load_init;           % Demand
            tmptmp_total_C = tmp_total_C; % Survived Index
            while redundancy_true

                % Compare the capacity and demand
                g_init = capacity-demand; % check failure

                % Determine whether fail or not
                if min(g_init) > 0
                    non_fail = non_fail+1;
                    redundancy_true = false;
                else
                    % Check cascading failure
                    if length(g_init)==1
                        redundancy_true = false;
                        fail = fail+1;
                    else
                                                
                        tmp_total_C2 = tmptmp_total_C; % Index of survived element

                        % Find the failed indx
                        tmp_idx = find(g_init == min(g_init));
                        if length(tmp_idx) > 1
                            tmp_idx = tmp_idx(1);
                        end
                        tmp_total_C2(tmp_idx) = []; % Erase the failed element

                        % New force for each element (demand)
                        demand_2 = ones([1,length(g_init) - 1])*force_mean/ ...
                                                                (length(g_init) - 1);

                        % Corresponding capacity
                        capacity_2 = capacity;
                        capacity_2(tmp_idx) = [];

                        % Change the capacity, demand, and index
                        capacity = capacity_2;
                        demand = demand_2;
                        tmptmp_total_C = tmp_total_C2;
                    end
                end

            end % end of while (check redundancy)

            %tmp_prob_history(kk) = fail/MCS_N;
        end % end of MCS for loop
    end % end of isempty    

    total_redu_prob = [total_redu_prob;[fail, non_fail]];
    disp('The current for loop ii is')
    disp(ii)    
end

redundancy_prob = total_redu_prob/MCS_N;


Redun_index = [];
for ii =1:length(redundancy_prob)
    tmp = redundancy_prob(ii,1);
    if -norminv(tmp)<-3
        Redun_index = [Redun_index; -3];
    else
        Redun_index = [Redun_index; -norminv(tmp)];
    end
        
end    

%% Plot (Figure 5)
% Plot scatter
color_code = [{'#0072BD'},	'#D95319', '#7E2F8E', '#000000', '#77AC30', 	'#A2142F'];
line_code = [{'-'},'--','-','-','-','--'];

close all
figure('Renderer', 'painters', 'Position', [10 10 700 700])

s1 = scatter(Redun_index, Relia_index,200, 'bo'); hold on;

xlim([-3,7]);ylim([-3,7]);
xticks([-3,-2,-1,0,1,2,3,4,5,6,7])
yticks([-3,-2,-1,0,1,2,3,4,5,6,7])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',30)
xlabel('Redundancy, \pi','FontSize',29); ylabel('Reliability, \beta','FontSize',29)
grid on


% Plot contour
struc_integrity = [10^-1, 10^-2, 10^-3, 10^-4, 10^-5];

for ii =1:length(struc_integrity)
    
    beta = -3:0.0001:-norminv(struc_integrity(ii))+1;
    gamma = -norminv(struc_integrity(ii)./normcdf(-beta));
    
    R2_value{ii} = [beta', gamma'];
end

for ii=4:4%length(struc_integrity)
    tmp = R2_value{ii};
    beta = tmp(:,1);
    gamma = tmp(:,2);
    tmp_color = sscanf(color_code{ii}(2:end),'%2x%2x%2x',[1 3])/255;

    eval(sprintf('p%d = plot(gamma, beta)', ii));
    eval(sprintf('p%d.LineWidth = 3', ii));
    eval(sprintf('p%d.Color = tmp_color', ii));
    eval(sprintf('p%d.LineStyle = line_code{ii}', ii));
    
    hold on
end
hold off

