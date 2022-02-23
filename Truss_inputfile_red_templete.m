clear probdata femodel analysisopt gfundata randomfield systems results output_filename


% Define the name and type of the file to store the output of the analysis:
% output_filename = 'outputfile_truss.txt';

random = 0;

%%% FINDME %%%

% scenario_red = {[]; [8 13 14 18 19]; [9 14 15 20 21]; ...
%     [10 15 16 22 23]; [11 16 17 24 25]; []; []; [2 13 14 18 19];...
%     [3 14 15 20 21]; [4 15 16 22 23]; [5 16 17 24 25]; []; ...
%     [1 2 7 8 14 18 19]; [2 3 8 9 13 14 15 18 19 20 21];[3 4 9 10 14 16 20 21 22 23]; ...
%     [4 5 10 11 15 17 22 23 24 25]; [5 6 11 12 24 25]; [2 8 13 14 19]; [2 8 13 14 18]; ...
%     [3 9 14 15 21]; [3 9 14 15 20]; [4 10 15 16 23]; [4 10 15 16 22]; [5 11 16 17 25]; ...
%     [5 11 16 17 24]};
% for ii = 1:length(scenario_red{scenario_rel})
%     if scenario_red{scenario_rel}(ii) > scenario_rel
%         scenario_red{scenario_rel}(ii) = scenario_red{scenario_rel}(ii)-1;
%     end
% end

probdata.marg(1,:) =  [ 1  -120*1e3       120*1e3*0.3    -120*1e3   0 0 0 0 0];
probdata.marg(2,:) =  [ 1  -140*1e3       140*1e3*0.3    -140*1e3   0 0 0 0 0];
probdata.marg(3,:) =  [ 1  -140*1e3       140*1e3*0.3    -140*1e3   0 0 0 0 0];
probdata.marg(4,:) =  [ 1  -120*1e3       120*1e3*0.3    -120*1e3   0 0 0 0 0];
probdata.marg(5,:) =  [ 1  -120*1e3       120*1e3*0.3    -120*1e3   0 0 0 0 0];

probdata.correlation = eye(size(probdata.marg,1)); 

probdata.parameter = distribution_parameter(probdata.marg);

analysisopt.ig_max     = 30;
analysisopt.il_max     = 5;
analysisopt.e1        = 0.001;
analysisopt.e2        = 0.001;
analysisopt.step_code = 0;
analysisopt.grad_flag = 'ffd';
analysisopt.sim_point = 'dspt';
analysisopt.stdv_sim  = 1;
analysisopt.num_sim   = 10000;
analysisopt.target_cov = 0.025;

femodel.ndf = 2;

femodel.node( 1,:) = [   0.0   0.0  ];
femodel.node( 2,:) = [   2.0   1.6  ];
femodel.node( 3,:) = [   4.0   2.2  ];
femodel.node( 4,:) = [   6.0   2.6  ];
femodel.node( 5,:) = [   8.0   2.2  ];
femodel.node( 6,:) = [  10.0   1.6  ];
femodel.node( 7,:) = [  12.0   0.0  ];
femodel.node( 8,:) = [   2.0   0.0  ];
femodel.node( 9,:) = [   4.0   0.0  ];
femodel.node(10,:) = [   6.0   0.0  ];
femodel.node(11,:) = [   8.0   0.0  ];
femodel.node(12,:) = [  10.0   0.0  ];


femodel.el(1,:) = [ 1  1 2 200e9  area(1) ];
femodel.el(2,:) = [ 1  2 3 200e9  area(2) ];
femodel.el(3,:) = [ 1  3 4 200e9  area(3) ];
femodel.el(4,:) = [ 1  4 5 200e9  area(4) ];
femodel.el(5,:) = [ 1  5 6 200e9  area(5) ];
femodel.el(6,:) = [ 1  6 7 200e9  area(6) ];
femodel.el(7,:) = [ 1  1 8 200e9  area(7) ];
femodel.el(8,:) = [ 1  8 9 200e9  area(8) ];
femodel.el(9,:) = [ 1  9 10 200e9  area(9) ];
femodel.el(10,:) =[ 1  10 11 200e9  area(10) ];
femodel.el(11,:) =[ 1  11 12 200e9  area(11) ];
femodel.el(12,:) =[ 1  12 7 200e9  area(12) ];
femodel.el(13,:) =[ 1  8 2 200e9  area(13) ];
femodel.el(14,:) =[ 1  9 3 200e9  area(14) ];
femodel.el(15,:) =[ 1  10 4 200e9  area(15) ];
femodel.el(16,:) =[ 1  11 5 200e9  area(16) ];
femodel.el(17,:) =[ 1  12 6 200e9  area(17) ];
femodel.el(18,:) =[ 1  8 3 200e9  area(18) ];
femodel.el(19,:) =[ 1  2 9 200e9  area(19) ];
femodel.el(20,:) =[ 1  9 4 200e9  area(20) ];
femodel.el(21,:) =[ 1  3 10 200e9  area(21) ];
femodel.el(22,:) =[ 1  4 11 200e9  area(22) ];
femodel.el(23,:) =[ 1  10 5 200e9  area(23) ];
femodel.el(24,:) =[ 1  5 12 200e9  area(24) ];
femodel.el(25,:) =[ 1  11 6 200e9  area(25) ];

if length(scenario_rel) == 1
if scenario_rel == 1
    femodel.el = femodel.el(2:end,:);
elseif scenario_rel == 25
    femodel.el = femodel.el(1:end-1,:);
else
    femodel.el = [femodel.el(1:scenario_rel-1,:);femodel.el(scenario_rel+1:end,:)];
end
elseif length(scenario_rel) == 2
    if scenario_rel(2) == 25
        femodel.el = femodel.el(1:end-1,:);
        if scenario_rel(1) == 1
            femodel.el = femodel.el(2:end,:);
        else
            femodel.el = [femodel.el(1:scenario_rel(1)-1,:);femodel.el(scenario_rel(1)+1:end,:)];
        end
    else
        femodel.el = [femodel.el(1:scenario_rel(2)-1,:);femodel.el(scenario_rel(2)+1:end,:)];
        if scenario_rel(1) == 1
            femodel.el = femodel.el(2:end,:);
        else
            femodel.el = [femodel.el(1:scenario_rel(1)-1,:);femodel.el(scenario_rel(1)+1:end,:)];
        end
    end
end

femodel.loading(1,:) = [ 8 (random)  2 1 ];
femodel.loading(2,:) = [ 9 (random)  2 1 ];
femodel.loading(3,:) = [ 10 (random)  2 1 ];
femodel.loading(4,:) = [ 11 (random)  2 1 ];
femodel.loading(5,:) = [ 12 (random)  2 1 ];

femodel.nodal_spring = 0;

femodel.fixed_dof(1,:) = [ 1  1 1 ];
femodel.fixed_dof(2,:) = [ 7  0 1 ];

for i = 1 : 5
   femodel.id(i,:) = [ i 1 i ];
end

gfundata(1).evaluator = 'FERUMlinearfecode';
gfundata(1).type = 'displacementlimit';
gfundata(1).parameter = 'yes';


gfundata(1).resp = [ 5 1 ];


gfundata(1).lim = 276e6;

randomfield.mesh = 0;

