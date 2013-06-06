function [max_lambda, predicted_list, corrected_list, combined_list, success, et] = cpf(casedata, participation, sigmaForLambda, sigmaForVoltage, verbose)
%CPF  Run continuation power flow (CPF) solver.
%   [INPUT PARAMETERS]
%   loadvarloc: load variation location(in external bus numbering). Single
%               bus supported so far.
%   sigmaForLambda: stepsize for lambda
%   sigmaForVoltage: stepsize for voltage
%   [OUTPUT PARAMETERS]
%   max_lambda: the lambda in p.u. w.r.t. baseMVA at (or near) the nose
%               point of PV curve
%   NOTE: the first column in return parameters 'predicted_list,
%   corrected_list, combined_list' is bus number; the last row is lambda.
%   created by Rui Bo on 2007/11/12

%   MATPOWER
%   $Id: cpf.m,v 1.7 2010/04/26 19:45:26 ray Exp $
%   by Rui Bo
%   and Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2010 by Power System Engineering Research Center (PSERC)
%   Copyright (c) 2009-2010 by Rui Bo
%
%   This file is part of MATPOWER.
%   See http://www.pserc.cornell.edu/matpower/ for more info.
%
%   MATPOWER is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation, either version 3 of the License,
%   or (at your option) any later version.
%
%   MATPOWER is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with MATPOWER. If not, see <http://www.gnu.org/licenses/>.
%
%   Additional permission under GNU GPL version 3 section 7
%
%   If you modify MATPOWER, or any covered work, to interface with
%   other modules (such as MATLAB code and MEX-files) available in a
%   MATLAB(R) or comparable environment containing parts covered
%   under other licensing terms, the licensors of MATPOWER grant
%   you additional permission to convey the resulting work.

%% define named indices into bus, gen, branch matrices

%escalate 'singular' to a matrix so we can use error handling to deal with
%it

[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, ...
    GEN_STATUS, PMAX, PMIN, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN] = idx_gen;

%% assign default parameters
if nargin < 3
    sigmaForLambda = 0.1;       % stepsize for lambda
    sigmaForVoltage = 0.025;    % stepsize for voltage
end

if nargin < 5, verbose = 0; end



%% options
max_iter = 1000;                 % depends on selection of stepsizes

%% ...we use PV curve slopes as the criteria for switching modes
% slopeThresh_Phase1 = 0.5;       % PV curve slope shreshold for voltage prediction-correction (with lambda increasing)
% slopeThresh_Phase2 = 0.3;       % PV curve slope shreshold for lambda prediction-correction
slopeThresh = [0.5, 0.3, Inf];

%% load the case & convert to internal bus numbering
[caseData.baseMVA, caseData.bus, caseData.gen, branch] = loadcase(casedata);
[numBuses, ~] = size(caseData.bus);

if nargin < 2 %no participation factors so keep the current load profile
	participation = caseData.bus(:,PD)./sum(caseData.bus(:,PD));
else
	
end

participation = participation(:);%(:) forces column vector

if length(participation) ~= numBuses, %improper number of participations given
	if length(participation) == 1 && participation > 0 && isinteger(participation),%assume bus number is specified instead
		participation = (1:numBuses)'==participation;
	else
		fprintf('\t[Info]\tParticipation Factors improperly specified.\n\t\t\tKeeping Current Loading Profile.\n');
		participation = caseData.bus(:,PD)./sum(caseData.bus(:,PD));
	end
end

[i2e, caseData.bus, caseData.gen, branch] = ext2int(caseData.bus, caseData.gen, branch);
e2i = sparse(max(i2e), 1);
e2i(i2e) = (1:size(caseData.bus, 1))';
%loadvarloc_i = e2i(loadvarloc);
participation_i = participation(e2i(i2e));

participation_i = participation_i ./ sum(participation_i); %normalize
caseData.participation = participation_i;

%% get bus index lists of each type of bus
[caseData.ref, caseData.pv, caseData.pq] = bustypes(caseData.bus, caseData.gen);

%% generator info
on = find(caseData.gen(:, GEN_STATUS) > 0);      %% which generators are on?
gbus = caseData.gen(on, GEN_BUS);                %% what buses are they at?

%% form Ybus matrix
[caseData.Ybus, ~, ~] = makeYbus(caseData.baseMVA, caseData.bus, branch);


%% initialize parameters
% set lambda to be increasing
flag_lambdaIncrease = true;  % flag indicating lambda is increasing or decreasing

%get all QP ratios
caseData.initQPratio = caseData.bus(:,QD)./caseData.bus(:,PD);
if any(isnan(caseData.initQPratio)), 
	if verbose > 1, fprintf('\t[Warning]:\tLoad real power at bus %d is 0. Q/P ratio will be fixed at 0.\n', find(isnan(initQPratio)));  end
	caseData.initQPratio(isnan(caseData.initQPratio)) = 0;
end


lambda0 = 0; 
lambda = lambda0;
Vm = ones(size(caseData.bus, 1), 1);          %% flat start
Va = caseData.bus(caseData.ref(1), VA) * Vm;
V  = Vm .* exp(sqrt(-1) * pi/180 * Va);
V(gbus) = caseData.gen(on, VG) ./ abs(V(gbus)).* V(gbus);

pointCnt = 0;

%% do voltage correction (ie, power flow) to get initial voltage profile
lambda_predicted = lambda;
V_predicted = V;
[V, lambda, success, ~] = cpf_correctVoltage(caseData.baseMVA, caseData.bus, caseData.gen, caseData.Ybus, V_predicted, lambda_predicted, caseData.initQPratio, caseData.participation);
%% record data
pointCnt = pointCnt + 1;

%%------------------------------------------------
% do cpf prediction-correction iterations
%%------------------------------------------------
t0 = clock;


states = {'voltage prediction-correction (lambda increasing)',...
		  'lambda prediction-correction (voltage decreasing)',...
		  'voltage prediction-correction (lambda decreasing)'};
	  
state = 1; %three states: lambda increase, voltage decrease, lambda increase

predicted_list = [];
corrected_list = [];

Vpr = [];
lambda_pr = [];
V_corr = [];
lambda_corr = [];

continuationBus = 0; %not used until after it is set
i = 0;



if verbose > 0  fprintf('Start Phase %d: %s.\n', state, states{state}); end

while state < length(states) && i < max_iter
	
	
	
%% --- Start Phase 1: voltage prediction-correction (lambda increasing)




    i = i + 1;
    
    % save good data
    V_saved = V;
    lambda_saved = lambda;
    
    %% do voltage prediction to find predicted point (predicting voltage)
    [V_predicted, lambda_predicted, J] = cpf_predict(caseData.Ybus, caseData.ref, caseData.pv, caseData.pq, V, lambda, sigmaForLambda, 1, caseData.initQPratio, caseData.participation, flag_lambdaIncrease);
    
    %% do voltage correction to find corrected point
	
	[V, lambda, success, ~] = cpf_correct(caseData, state, V_predicted, lambda_predicted, continuationBus);
	

    %% calculate slope (dP/dLambda) at current point
	[slope, continuationBus] = max(abs(V-V_saved)  ./ (lambda-lambda_saved)); %calculate maximum slope at current point.
 
	
	
    %% instead of using condition number as criteria for switching between
    %% modes...
    %%    if rcond(J) <= condNumThresh_Phase1 | success == false % Jacobian matrix is ill-conditioned, or correction step fails
    %% ...we use PV curve slopes as the criteria for switching modes
    if (state == 1 && abs(slope) >= slopeThresh(state)) || (state == 2 && abs(slope) <= slopeThresh(state)) || success == false || lambda < 0 % Approaching nose area of PV curve, or correction step fails
        % restore good data

		lambda = lambda_saved;
        
		if lambda < 0 % lambda is less than 0, then stops CPF simulation
			if verbose > 0
				fprintf('\t[Info]:\tlambda is less than 0.\n\t\t\tCPF finished.\n');
			end
			state = state+1; continue;
		end
		if verbose > 0
			if ~success, 
				if ~isempty(strfind(lastwarn, 'singular')), fprintf('\t[Info]:\tMatrix is singular. Aborting Correction.\n'); break;
				else fprintf('\t[Info]:\tLambda correction fails.\n'); end
			else
				switch state,
					case 1,	fprintf('\t[Info]:\tApproaching nose area of PV curve.\n');
					case 2, fprintf('\t[Info]:\tLeaving nose area of PV curve.\n');
					case 3, fprintf('\t[Info]:\tFinished.\n');
				end
			end
		end
         
		if state==3, flag_lambdaIncrease = false; end
		if verbose > 0,	fprintf('\t[Info]:\t%d data points contained in phase %d.\n', i-1, state); end
		state = state+1;
		i = 0;
		if verbose > 0,	fprintf('Switch to Phase %d: %s.\n', state, states{state});	end
		
		continue;
    else
        if verbose > 2
            fprintf('\nVm_predicted\tVm_corrected\n');
            [[abs(V_predicted);lambda_predicted] [abs(V);lambda]]
        end

        %% record data
	
		
        pointCnt = pointCnt + 1;
		Vpr = [Vpr, V_predicted];
		lambda_pr = [lambda_pr, lambda_predicted];
		V_corr = [V_corr, V];
		lambda_corr = [lambda_corr, lambda];
		
        predicted_list(:, pointCnt-1) = [V_predicted;lambda_predicted];
        corrected_list(:, pointCnt) = [V;lambda];
    end
end


et = etime(clock, t0);

if i == max_iter,
	fprintf('\t[Info] Max iterations hit.\n');
end

if  ~exist('corrected_list','var') | ~exist('predicted_list','var'),

	max_lambda = NaN;
	predicted_list = [];
	corrected_list = [];
	combined_list = [];
	success==false;
end

if ~exist('max_lambda', 'var'),
	max_lambda = NaN;
	success==false;
end

if ~isempty(predicted_list) && ~isempty(corrected_list),
	%% combine the prediction and correction data in the sequence of appearance
	% NOTE: number of prediction data is one less than that of correction data
	predictedCnt = size(predicted_list, 2);
	combined_list(:, 1) = corrected_list(:, 1);
	for i = 1:predictedCnt
		combined_list(:, 2*i)     = predicted_list(:, i);
		combined_list(:, 2*i+1)   = corrected_list(:, i+1);
	end

	%% convert back to original bus numbering & print results
	[caseData.bus, caseData.gen, branch] = int2ext(i2e, caseData.bus, caseData.gen, branch);

	%% add bus number as the first column to the prediction, correction, and combined data list
	nb          = size(caseData.bus, 1);
	max_lambda  = max(corrected_list(nb+1, :));
	predicted_list = [[caseData.bus(:, BUS_I);0] predicted_list];
	corrected_list = [[caseData.bus(:, BUS_I);0] corrected_list];
	combined_list  = [[caseData.bus(:, BUS_I);0] combined_list];
else
	combined_list = [];
end

if verbose > 1
    Vm_corrected = abs(corrected_list);
    Vm_predicted = abs(predicted_list);
    Vm_combined  = abs(combined_list);
    Vm_corrected
    Vm_predicted
    Vm_combined
    pointCnt_Phase1
    pointCnt_Phase2
    pointCnt_Phase3
    pointCnt
end

if nargout == 1,
	results.max_lambda = max_lambda;
	results.V_pr = Vpr;
	results.lambda_pr = lambda_pr;
	results.V_corr = V_corr;
	results.lambda_corr = lambda_corr;
	results.success = success;
	results.et = et;
	
	max_lambda = results; %return a single struct
end

%% Changelog - Anton Lodder - 2013.3.27
% I implemented participation factor loading to allow all buses to
% participate in load increase as a function of lambda.
%
% * Participation factors should be given as a vector, one value for each
%   bus.
% * if only one value is given, it is assumed that the value is a bus
%   number rather than a  participation factor, and all other buses get a
%   participation factor of zero (point two)
% * any buses with zero participation factor will remain at their initial
%   load level
% * from the previous two bullets: backwards compatibility is maintained
%   while allowing increased functionality
% * if participation is not a valid bus number (eg float, negative number),
%	maintains given bus loading profile.
% * if no participation factors are given, maintain given bus loading
%   profile.
