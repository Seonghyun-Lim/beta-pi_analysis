

clear;  clc;

%% output vector setting

all_scenarios = [1:25]';    % number of failure scenarios

beta = zeros(25,2);   P_fb = zeros(25,2);

pi_r = zeros(25,2);
beta_sing = zeros(length(all_scenarios),25);  % beta^(i)_j in Eq.(9)
P_fd = zeros(25,2);
alpha_sing = zeros(5,25,length(all_scenarios));


for idx_rein = 1:3    % 1(No reinforcement), 2(Reinforcement 1), 3(Reinforcement 2)
    % Geometry Setting   % Area of members which determine structure resilience
A = 32e-4;  B = 24e-4;  C = 5e-4;  D = 7e-4;  % four types of cross-section (m^2)
ar = [A A A A A A B B B B B B C C C C C D D D D D D D D]; % area of cross-section
if idx_rein ==2
    ar = ar + [0 0 0 0 0 0 0 0 0 0 0 0 3 0 0 0 3 0 0 0 0 0 0 0 0]*1e-4;  % reinforcement to elements #13 and #17
elseif idx_rein ==3
    ar = ar + [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 0 0 0 0 0 3]*1e-4;  % reinforcement to elements #19 and #25
end

%% Estimate reliability index
for scenario = 1:25
    % generating ferum input file
    fid1 = fopen('Truss_inputfile_rel_templete.m','r');
    lines_ = textscan(fid1,'%s','delimiter','\n');
    lines=lines_{1};
    fclose(fid1);
    for il=1:numel(lines)
        if isempty(regexp(lines{il},'%%% FINDME %%%'))==0
            FINDME=il;
        end
    end
    copyfile Truss_inputfile_rel_templete.m Truss_inputfile_rel.m;
    fid2 = fopen('Truss_inputfile_rel.m', 'w');
    for nl = 1:numel(lines)
        fprintf(fid2,'%s \n', lines{nl});
        if nl == FINDME
            fprintf(fid2,'gfundata(1).scenario = [%i];\n',scenario);
            fprintf(fid2,'area = [%i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i];'...
                ,ar(1),ar(2),ar(3),ar(4),ar(5),ar(6),ar(7),ar(8),ar(9),ar(10),ar(11),ar(12),...
                ar(13),ar(14),ar(15),ar(16),ar(17),ar(18),ar(19),ar(20),ar(21),ar(22),ar(23),ar(24),ar(25));
        end
    end
    fclose(fid2);
    run('Truss_inputfile_rel.m')
    
    ferum;
    beta(scenario,idx_rein) = real(formresults.beta1);
    P_fb(scenario,idx_rein) = formresults.pf1;
end

%% Estimate redundancy index
epsilon = 1e-10;
fid1 = fopen('Truss_inputfile_red_templete.m','r'); 
lines_ = textscan(fid1,'%s','delimiter','\n');
lines = lines_{1};
fclose(fid1);
for il = 1:numel(lines)
    if isempty(regexp(lines{il},'%%% FINDME %%%'))==0
        FINDME = il;
    end
end

for ii = [2:5 8:11 13:25]   % Every failure scenario considered in redundancy phase
    scenario = all_scenarios(ii,:);
    num_comp_red = 24;  % eliminating failed components
    for jj = 1:length(ar)  % remaining components
        % generating ferum input file
        copyfile Truss_inputfile_red_templete.m Truss_inputfile_red.m;
        fid2 = fopen('Truss_inputfile_red.m', 'w');
        kk = jj;   % kk is real no. of component
        if jj >= scenario
            jj = jj-1;  % jj is shrinked no. of component
        end
        for nl = 1:numel(lines)
            fprintf(fid2,'%s \n', lines{nl});
            if nl == FINDME
                fprintf(fid2,'scenario_rel = %i;\n',scenario);
                if kk == 1
                    fprintf(fid2,'gfundata(1).scenario = [%i];\n',kk);
                else
                    fprintf(fid2,'gfundata(1).scenario = [%i];\n',jj);
                end
                fprintf(fid2,'area = [%i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i];'...
                    ,ar(1),ar(2),ar(3),ar(4),ar(5),ar(6),ar(7),ar(8),ar(9),ar(10),ar(11),ar(12),...
                    ar(13),ar(14),ar(15),ar(16),ar(17),ar(18),ar(19),ar(20),ar(21),ar(22),ar(23),ar(24),ar(25));
            end
        end
        fclose(fid2);
        run('Truss_inputfile_red.m')
        
        ferum;
        
        % beta_sing: beta^(i)_j in Eq.(9)
        if kk == scenario
            beta_sing(ii,kk) = NaN;    
        else
            beta_sing(ii,kk) = real(formresults.beta1);
        end
        if (kk == scenario) || (sum(ii==[1 6 7 12])==1)
            alpha_sing(:,kk,ii) = NaN;
        else
            alpha_sing(:,kk,ii) = formresults.alpha;
        end
    end
    
    % correlation coefficient matrix calculation 
    corrmat = zeros(24,24);     % correlation coefficient matrix in Eq.(9)
    if ii == 25
        alpha_corr = alpha_sing(:,1:24,ii);
    else
        alpha_corr = alpha_sing(:,[1:ii-1,ii+1:end],ii);
    end
    
    for mm = 1:24
        for nn = 1:24
            if mm == nn
                corrmat(mm,nn) = 1;
            else
                corrmat(mm,nn) = alpha_corr(:,mm)'*alpha_corr(:,nn);
            end
        end
    end
    try chol(corrmat)
        disp('Matrix is symmetric positive definite.')
    catch ME
        disp('Matrix is not symmetric positive definite')
        corrmat = 1/(1+epsilon)*(corrmat+diag(repmat(epsilon,1,24)));
    end
    
    P_system = 1 - mvncdf(rmmissing(beta_sing(ii,:)),zeros(1,24),corrmat);  % Eq.(9)
    pi_r(ii,idx_rein) = -norminv(P_system);   % Eq.(2)
end
clc;
end

%% Plot (Figure 11)
colorcode_yellow = [235 149 3]/255;

pi_r([1,6,7,12],:) = [-3,-3,-3; -3,-3,-3; -3,-3,-3; -3,-3,-3]; 
beta(beta == -Inf) = -3;   % left bound of beta 
pi_r(pi_r == -Inf) = -3;   % left bound of pi

[xct1, xct2] = meshgrid(-3:0.05:5, -3:0.05:5);
sys_integ = normcdf(-xct1).*normcdf(-xct2);

figure('Renderer', 'painters', 'Position', [10 10 880 1000])   % Figure 11
scatter(pi_r(:,1),beta(:,1),150,'bo','LineWidth', 1.5)
hold on
scatter(pi_r(:,2),beta(:,2),150,'rx','LineWidth', 1.5)
scatter(pi_r(:,3),beta(:,3),150,colorcode_yellow,'+','LineWidth', 1.5)
plot([100 100],[100 101],'LineStyle',':','LineWidth', 3,'color','k')
plot([100 100],[100 101],'LineStyle','-','LineWidth', 3,'color','k')
ar = get(gca,'XTickLabel');
set(gca,'XTickLabel',ar,'FontName','Times','fontsize',30)

contour(xct1, xct2, sys_integ,'levellist',1e-4,'LineStyle',':','LineWidth', 3,'color','k');
xlim([-3 7]); ylim([-3 9])
xticks(-3:2:7); yticks(-3:2:9)
ar = get(gca,'XTickLabel');
set(gca,'XTickLabel',ar,'FontName','Times','fontsize',30)
xlabel('Redundancy, \pi','FontSize',29); ylabel('Reliability, \beta','FontSize',29)
daspect([1 1 1])
[h,icons,plots,legend_text] = legend('No reinforcement','Reinforcement 1','Reinforcement 2','FontSize',20)
icons(6).Children.MarkerSize = 12;
icons(7).Children.MarkerSize = 12;
icons(8).Children.MarkerSize = 12;
icons(6).Children.LineWidth = 1.5;
icons(7).Children.LineWidth = 1.5;
icons(8).Children.LineWidth = 1.5;
grid on
hold off


