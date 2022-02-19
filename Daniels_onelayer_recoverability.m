
clear; clc; close all;

%% Set parameters
rng(117)
num_bar = 15;

force_mean = (num_bar)*200; % External force kN
Statistical_data = ones(num_bar); % Area information, mm2
cap_cov = 0.35;
MCS_N = 1e6;

%% Estimate reliability index
% Daniels system consists of 3 types of member: 'S(320Mpa),'M(400Mpa)','L(500Mpa)'
cap_meanS = 320;
cap_meanM = 400;
cap_meanL = 500;
sigS = cap_meanS*cap_cov;
sigM = cap_meanM*cap_cov;
sigL = cap_meanL*cap_cov;

% reliability index for a single component failure
betaS = (cap_meanS - force_mean/num_bar)/sigS;
betaM = (cap_meanM - force_mean/num_bar)/sigM;
betaL = (cap_meanL - force_mean/num_bar)/sigL;
% failure probability for a single component failure
pfS = normcdf(-betaS);
pfM = normcdf(-betaM);
pfL = normcdf(-betaL);

% double-component-failure scenario
scenario_name_2 = {[num2str(cap_meanS),', ',num2str(cap_meanS)] [num2str(cap_meanS),', ',num2str(cap_meanM)]...
    [num2str(cap_meanS),', ',num2str(cap_meanL)] [num2str(cap_meanM),', ',num2str(cap_meanM)]...
    [num2str(cap_meanM),', ',num2str(cap_meanL)] [num2str(cap_meanL),', ',num2str(cap_meanL)]};



relpSS = pfS^2*(1-pfS)^3*(1-pfM)^5*(1-pfL)^5;       % failure of two 'S'
relpSM = pfS*pfM*(1-pfS)^4*(1-pfM)^4*(1-pfL)^5;     % failure of one 'S' and one 'M'
relpSL = pfS*pfL*(1-pfS)^4*(1-pfM)^5*(1-pfL)^4;     % failure of one 'S' and one 'L'
relpMM = pfM^2*(1-pfM)^3*(1-pfL)^5*(1-pfS)^5;       % failure of two 'M'
relpML = pfM*pfL*(1-pfS)^5*(1-pfM)^4*(1-pfL)^4;     % failure of one 'M' and one 'L'
relpLL = pfL^2*(1-pfS)^5*(1-pfM)^5*(1-pfL)^3;       % failure of two 'L'


relp_vec = [relpSS, relpSM, relpSL, relpMM, relpML, relpLL]';   % failure probability of double-component-failure scenarios
relia_index_2 = -norminv(relp_vec);   % reliability indices


%% Estimate redundancy index

capacity_remain = [repmat(cap_meanS,1,3), repmat(cap_meanM,1,5) repmat(cap_meanL,1,5);...
    repmat(cap_meanS,1,4), repmat(cap_meanM,1,4) repmat(cap_meanL,1,5);...
    repmat(cap_meanS,1,4), repmat(cap_meanM,1,5) repmat(cap_meanL,1,4);...
    repmat(cap_meanS,1,5), repmat(cap_meanM,1,3) repmat(cap_meanL,1,5);...
    repmat(cap_meanS,1,5), repmat(cap_meanM,1,4) repmat(cap_meanL,1,4);...
    repmat(cap_meanS,1,5), repmat(cap_meanM,1,5) repmat(cap_meanL,1,3)];

total_redu_prob = [];

for scenario = 1:6
    
    ii = 2;
    % Failure cases where some elements are survived
    % Load for each survived element
    load_init = ones([1,num_bar-ii])*force_mean/(num_bar-ii);
    tmp_total_C = 1:1:length(Statistical_data)-ii;
    
    % MCS is carried out
    non_fail = 0; fail = 0; %tmp_prob_history = zeros(MCS_N,1);
    for kk = 1 : MCS_N
        
        % Randomly generate capacities
        capacity = normrnd(capacity_remain(scenario,:), capacity_remain(scenario,:)*cap_cov);
        
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
    % end of isempty

total_redu_prob = [total_redu_prob;[fail, non_fail]];
disp('The current for loop ii is')
disp(ii)
end

redundancy_prob = total_redu_prob/MCS_N;


redun_index_2 = [];
for ii =1:length(redundancy_prob)
    tmp = redundancy_prob(ii,1);
    if -norminv(tmp)<-3
        redun_index_2 = [redun_index_2; -3];
    else
        redun_index_2 = [redun_index_2; -norminv(tmp)];
    end
    
end

%% Incorporating recoverability
% assume that the recovery cost is the sum of the exponential functions of the yield stresses of the failed members
recovSS = exp(cap_meanS/100)+exp(cap_meanS/100);
recovSM = exp(cap_meanS/100)+exp(cap_meanM/100);
recovSL = exp(cap_meanS/100)+exp(cap_meanL/100);
recovMM = exp(cap_meanM/100)+exp(cap_meanM/100);
recovML = exp(cap_meanM/100)+exp(cap_meanL/100);
recovLL = exp(cap_meanL/100)+exp(cap_meanL/100);

recover_index_2 = [recovSS recovSM recovSL recovMM recovML recovLL];


%% Plot (Figure 6)
% Plot scatter
close all
figure('Renderer', 'painters', 'Position', [10 10 1000 800])

s1 = scatter(redun_index_2, relia_index_2,300, recover_index_2,'filled'); hold on;
text(redun_index_2([1:3,5,6])+0.05,relia_index_2([1:3,5,6]),scenario_name_2([1:3,5,6]),'FontName','Times','FontSize',18)
text(redun_index_2(4)-0.20,relia_index_2(4),scenario_name_2(4),'FontName','Times','FontSize',18)
colormap jet
c = colorbar('Ticks',[70 275],...
         'TickLabels',{'High','Low'},'Direction','reverse');
c.Label.String = 'Recoverability';
c.Label.FontSize = 20;
c.FontSize = 20;

xlim([0 1.4]); ylim([2.2 3.6])
xticks([0 0.5 1]); yticks([2.5 3 3.5])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',30)
xlabel('Redundancy, \pi','FontSize',29); ylabel('Reliability, \beta','FontSize',29)
grid on


% Plot contour
struc_integrity = 10^-3.5;
beta = -3:0.0001:-norminv(struc_integrity)+1;
gamma = -norminv(struc_integrity./normcdf(-beta));

color_code = [{'#0072BD'},	'#D95319', '#7E2F8E', '#000000', '#77AC30', 	'#A2142F'];
line_code = [{'-'},'--','-','-','-','--'];
for ii=4
    tmp_color = sscanf(color_code{ii}(2:end),'%2x%2x%2x',[1 3])/255;
    
    eval(sprintf('p%d = plot(gamma, beta)', ii));
    eval(sprintf('p%d.LineWidth = 3', ii));
    eval(sprintf('p%d.Color = tmp_color', ii));
    eval(sprintf('p%d.LineStyle = line_code{ii}', ii));
    
    hold on
end
