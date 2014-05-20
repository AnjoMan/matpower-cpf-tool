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

% [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
%     VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus();
% 
% [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
%     RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST] = idx_brch();
% [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, ...
%     GEN_STATUS, PMAX, PMIN, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN] = idx_gen();

if nargin == 0 % run test suite
%     which('cpf.m')
    [path, name, ext] = fileparts(which(sprintf('%s.m',mfilename)));
    runFile = fullfile(path, 'test', sprintf('test_%s.m', name));
    
    fprintf('Running unit-test for <%s%s>  \n=============================\n\n',name,ext);
    run(runFile)
%     pwd
    
%     fullfile(path , 'test', 'test_cpf.m')
%     run(fullfile('test', 'test_cpf.m'));
    
    return;
end

BUS_I = 1;
PD = 3;
QD = 4;
VA = 9;

GEN_BUS = 1;
VG = 6;
GEN_STATUS = 8;



lastwarn('No Warning');
%% assign default parameters
if nargin < 3
    sigmaForLambda = 0.1;       % stepsize for lambda
    sigmaForVoltage = 0.025;    % stepsize for voltage
end

if nargin < 5, verbose = 0; end

if verbose, fprintf('CPF\n'); figure; end

mError = MException('CPF:cpf', 'cpf_error');
%% options
max_iter = 1000;                 % depends on selection of stepsizes

%% ...we use PV curve slopes as the criteria for switching modes
slopeThresh_Phase1 = 0.5;       % PV curve slope shreshold for voltage prediction-correction (with lambda increasing)
slopeThresh_Phase2 = 0.3;       % PV curve slope shreshold for lambda prediction-correction


%% load the case & convert to internal bus numbering
[baseMVA, busE, genE, branchE] = loadcase(casedata);
[numBuses, ~] = size(busE);

shouldIPlotEverything = false;
whatBusShouldIPlot = min(size(busE, 1), 11);

shouldIPlotEverything = true;

if nargin < 2 %no participation factors so keep the current load profile
	participation = busE(:,PD)./sum(busE(:,PD));
else
	
	participation = participation(:);%(:) forces column vector

	if length(participation) ~= numBuses, %improper number of participations given
		if length(participation) == 1 && participation > 0 && isinteger(participation),%assume bus number is specified instead
			participation = (1:numBuses)'==participation;
		else
			if verbose, fprintf('\t[Info]\tParticipation Factors improperly specified.\n\t\t\tKeeping Current Loading Profile.\n'); end
			participation = busE(:,PD)./sum(busE(:,PD));
		end
	end
end

[i2e, bus, gen, branch] = ext2int(busE, genE, branchE);
e2i = sparse(max(i2e), 1);
e2i(i2e) = (1:size(bus, 1))';
participation_i = participation(e2i(i2e));

participation_i = participation_i ./ sum(participation_i); %normalize


	
	
%% get bus index lists of each type of bus
[ref, pv, pq] = bustypes(bus, gen);



if any(isnan(participation_i)), %could happen if no busses had loads
	participation_i = zeros(length(participation_i), 1);
	participation_i(pq) = 1/numel(participation_i(pq));
end

%% generator info
on = find(gen(:, GEN_STATUS) > 0);      %% which generators are on?
gbus = gen(on, GEN_BUS);                %% what buses are they at?

%% form Ybus matrix

[Ybus, ~, ~] = makeYbus(baseMVA, bus, branch);
if det(Ybus) == 0
    addCause(mError, MException('MATPOWER:makeYBus', 'Ybus is singular'));
end

%% initialize parameters
% set lambda to be increasing
flag_lambdaIncrease = true;  % flag indicating lambda is increasing or decreasing

%get all QP ratios
initQPratio = bus(:,QD)./bus(:,PD);
if any(isnan(initQPratio)), 
	if verbose > 1, fprintf('\t[Warning]:\tLoad real power at bus %d is 0. Q/P ratio will be fixed at 0.\n', find(isnan(initQPratio)));  end
	initQPratio(isnan(initQPratio)) = 0;
end


lambda0 = 0; %changed from 2 -> 0 
lambda = lambda0;
Vm = ones(size(bus, 1), 1);          %% flat start
Va = bus(ref(1), VA) * Vm;
V  = Vm .* exp(1i* pi/180 * Va);
V(gbus) = gen(on, VG) ./ abs(V(gbus)).* V(gbus);

pointCnt = 0;

%% do voltage correction (ie, power flow) to get initial voltage profile
lambda_predicted = lambda;
V_predicted = V;
[V, lambda, success, iterNum] = cpf_correctVoltage(baseMVA, bus, gen, Ybus, V_predicted, lambda_predicted, initQPratio, participation_i);
%% record data
if any(isnan(V))
    err = MException('CPF:correctVoltageError', 'Generating initial voltage profile');
    errCause = MException('CPF:correctVoltageError', ['NaN bus voltage at ', mat2str(i2e(isnan(V)))]);
    err = addCause(err,errCause);
    
    throw(err);
end

pointCnt = pointCnt + 1;


%%------------------------------------------------
% do cpf prediction-correction iterations
%%------------------------------------------------
t0 = clock;
%% --- Start Phase 1: voltage prediction-correction (lambda increasing)
if verbose > 0
    fprintf('Start Phase 1: voltage prediction-correction (lambda increasing).\n');
end


Vpr = [];
lambda_pr = [];
V_corr = [];
lambda_corr = [];
correctionIters = []; 
slopes = [];





pointCnt = pointCnt + 1;
Vpr = [Vpr, V_predicted];
lambda_pr = [lambda_pr, lambda_predicted];
V_corr = [V_corr, V];
lambda_corr = [lambda_corr, lambda];
correctionIters = [correctionIters, iterNum];




% stepSize = 1

%parametrize step size for Phase 1
minStepSize = 0.01;
maxStepSize = 2;


stepSize = min(max(sigmaForLambda, minStepSize),maxStepSize);

finished = false;
continuationBus = pq(1);

i = 0; j=0; k=0; %initialize counters for each phase to zero
while i < max_iter
    %% update iteration counter
    i = i + 1;
    
    % save good data
    V_saved = V;
    lambda_saved = lambda;
    
	if i<4, 
		%% do voltage prediction to find predicted point (predicting voltage)
		[V_predicted, lambda_predicted, J] = cpf_predict(Ybus, ref, pv, pq, V, lambda, stepSize, 1, initQPratio, participation_i, flag_lambdaIncrease);
	else
		[V_predicted, lambda_predicted] = cpf_predict_voltage(V_corr, lambda_corr, lambda, stepSize, ref, pv, pq);
    end
    
           
    %% check prediction to make sure step is not too big so as to cause non-convergence
    
    error_predicted = max(abs(V_predicted- V));
    if ~success && error_predicted > 0.1
        newStepSize = 0.8*stepSize; %cut down the step size to reduce the prediction error
        if newStepSize > minStepSize,
            if verbose, fprintf('\t\tpredicted voltage change: %.2f. Step Size reduced from %.3f to %.3f\n', error_predicted, stepSize, newStepSize); end
            stepSize = newStepSize;
            i = i-1;            
            continue;
        end       
    end
    
    %% do voltage correction to find corrected point
    [V, lambda, success, iterNum] = cpf_correctVoltage(baseMVA, bus, gen, Ybus, V_predicted, lambda_predicted, initQPratio, participation_i);
	
    % if voltage correction fails, reduce step size and try again
    if success == false  && stepSize > minStepSize,		
		stepSize = stepSize * 0.3;		
		if verbose, fprintf('\t\tdidnt solve; changed stepsize to: %f\n', stepSize); end
		i = i-1;
		V = V_saved;
		lambda = lambda_saved;
		continue;
    end

    %% calculate slope (dP/dLambda) at current point
	[slope, ~] = max(abs(V-V_saved)  ./ (lambda-lambda_saved)); %calculate maximum slope at current point.
	

    
    %% if successful, check max error and adjust step sized to get a better error (meaning, balanced between a bigger step size and making sure we still converge)
	error = abs(V-V_predicted)./abs(V);    
	if abs(log(mean(error)/0.0001)) > 1 && mean(error)>0,
        newStepSize = stepSize - 0.1*log(mean(error)/0.0001); %adjust step size
        newStepSize = max( min(newStepSize,maxStepSize),minStepSize); %clamp step size
        
		if verbose, fprintf('\t\tmean prediction error: %.15f. changed stepsize from %.2f to %.2f\n', mean(error), stepSize, newStepSize); end
        stepSize = newStepSize;
	end
	

    
    if success % if correction converged we can save the point and do plotting/output in verbose mode
        if verbose > 2
            fprintf('\nVm_predicted\tVm_corrected\n');
            [[abs(V_predicted);lambda_predicted] [abs(V);lambda]]
        end

        %% record data
        logStepResults()
        pointCnt = pointCnt + 1;
    end
    
    
    	
    % instead of using condition number as criteria for switching between
    % modes...
    %    if rcond(J) <= condNumThresh_Phase1 | success == false % Jacobian matrix is ill-conditioned, or correction step fails
    % ...we use PV curve slopes as the criteria for switching modes:    
    if abs(slope) >= slopeThresh_Phase1 || success == false % Approaching nose area of PV curve, or correction step fails
        
        % restore good data point if convergence failed
        if success == false
            V = V_saved;
            lambda = lambda_saved;
        end 

		if verbose > 0
			if ~success, 
				if ~isempty(strfind(lastwarn, 'singular')), 
					fprintf('\t[Info]:\tMatrix is singular. Aborting Correction.\n'); 
					lastwarn('No error');
					break;
				else
					fprintf('\t[Info]:\tLambda correction fails.\n'); 
				end
			else
				fprintf('\t[Info]:\tApproaching nose area of PV curve.\n');
			end
		end
        break;    
    end
end


% fprintf('Average prediction error for voltage: %f\n', mean(mean( abs( V_corr - Vpr)./abs(V_corr))));
% fprintf('Avg num of iterations: %f\n', mean(correctionIters));
pointCnt_Phase1 = pointCnt; % collect number of points obtained at this phase
if verbose > 0
    fprintf('\t[Info]:\t%d data points contained in phase 1.\n', pointCnt_Phase1);
end








%% --- Switch to Phase 2: lambda prediction-correction (voltage decreasing)
if verbose > 0
    fprintf('Switch to Phase 2: lambda prediction-correction (voltage decreasing).\n');
end

correctionIters2 = [];
predictionErrors2 = [];


maxStepSize = 0.1;
minStepSize = 0.0001;
% maxStepSize_voltage = 0.1;
% minStepSize_voltage = 0.0001;
startSlope = slope;
stepSize = sigmaForVoltage;
j = 0;
while j < max_iter && ~finished
    %% update iteration counter
    j = j + 1;

    % save good data
    V_saved = V;
    lambda_saved = lambda;
    
    
    %% do lambda prediction to find predicted point (predicting lambda)
	if j<4
        [V_predicted, lambda_predicted, J] = cpf_predict(Ybus, ref, pv, pq, V, lambda, stepSize, [2, continuationBus], initQPratio, participation_i,flag_lambdaIncrease);
    else
        [V_predicted, lambda_predicted] = cpf_predict_lambda(V_corr, lambda_corr, lambda, stepSize, continuationBus, ref, pv, pq);
    end
   
    
    
    %% do lambda correction to find corrected point
    Vm_assigned = abs(V_predicted);
    
  
	
	[V, lambda, success, iterNum] = cpf_correctLambda(baseMVA, bus, gen, Ybus, Vm_assigned, V_predicted, lambda_predicted, initQPratio, participation_i, ref, pv, pq, continuationBus);
	
    mean_step = mean( abs(V_predicted-V_saved));    
	prediction_error = mean(abs(V-V_predicted)./abs(V));
    error_order = log(prediction_error/0.000001);
    
    if ( mean_step > 0.00001 && ~success) % if we jumped too far and correction step didn't converge
        newStepSize = stepSize * 0.4;
        if newStepSize > minStepSize, %if we are not below min step-size threshold go back and try again with new stepSize
            if verbose, fprintf('\t\tpredicted voltage change: %f. Step Size reduced from %.5f to %.5f\n', mean_step, stepSize, newStepSize); end
            stepSize = newStepSize;
            V = V_saved;
            lambda= lambda_saved;
            j =  j-1;
            continue;
        end

    elseif abs(error_order) > 2 && prediction_error>0, %if we havent just dropped the stepSize, consider changing it to get a better error outcome
%         newStepSize = stepSize - 0.1*log(prediction_error/0.000001); %adjust step size
% 	proposedStepSize = min( max(stepSize - 0.03*log( mean(prediction_error)/0.0001), minStepSize), maxStepSize);

        newStepSize = stepSize * 1.1;
        newStepSize = max( min(newStepSize,maxStepSize),minStepSize); %clamp step size
        
		if verbose, fprintf('\t\tmean prediction error: %.15f. changed stepsize from %.5f to %.5f\n', prediction_error, stepSize, newStepSize); end
        stepSize = newStepSize;
	end
	
	
	prediction_error = mean(abs(V-V_predicted)./abs(V));
	

    %% calculate slope (dP/dLambda) at current point
    mSlopes = abs(V-V_saved)./(lambda-lambda_saved);
    [~,continuationBusPQ] = max(abs(mSlopes(pq))); %limit to only PQ busses
    continuationBus = pq(continuationBusPQ);
    slope = mSlopes(continuationBus);
	slopes = [slopes slope];%log the slope

    
    if success %if correction step converged, log values and do verbosity
        if verbose > 2
            fprintf('\nVm_predicted\tVm_corrected\n');
            [[abs(V_predicted);lambda_predicted] [abs(V);lambda]]
        end

        %% record data
		logStepResults()
        pointCnt = pointCnt + 1;
    end
    
    
    
    %% instead of using condition number as criteria for switching between
    %% modes...
    %%    if rcond(J) >= condNumThresh_Phase2 | success == false % Jacobian matrix is good-conditioned, or correction step fails
    %% ...we use PV curve slopes as the criteria for switching modes:
    if ~success || (slope < 0 && slope > -slopeThresh_Phase2),

        % restore good data
        V = V_saved;
        lambda = lambda_saved;

        %% ---change to voltage prediction-correction (lambda decreasing)
        if verbose > 0
			if ~success, 
				if ~isempty(strfind(lastwarn, 'singular'))
					fprintf('\t[Info]:\tMatrix is singular. Aborting Correction.\n');
					lastwarn('No error');
					break;
				else
					fprintf('\t[Info]:\tLambda correction fails.\n');
				end
			else
				fprintf('\t[Info]:\tLeaving nose area of PV curve.\n');
			end
        end
        break;   
    end
    
%     if verbose, fprintf('lambda: %.3f,     slope: %.4f     error: %e    error_order: %f     stepSize: %.15f\n', lambda, slope, prediction_error, error_order,stepSize); end
end




% fprintf('Average prediction error for voltage: %f\n', mean(predictionErrors2));
% fprintf('Avg num of iterations: %f\n', mean(correctionIters2));

pointCnt_Phase2 = pointCnt - pointCnt_Phase1; % collect number of points obtained at this phase
if verbose > 0
    fprintf('\t[Info]:\t%d data points contained in phase 2.\n', pointCnt_Phase2);
end









%% --- Switch to Phase 3: voltage prediction-correction (lambda decreasing)
if verbose > 0
    fprintf('Switch to Phase 3: voltage prediction-correction (lambda decreasing).\n');
end
% set lambda to be decreasing
flag_lambdaIncrease = false; 


%set step size for Phase 3
maxIncrease = 0.3;
minStepSize = 0.01;
maxStepSize = 2;

k = 0;
while k < max_iter && ~finished
    %% update iteration counter
    k = k + 1;
    
    
    %% store V and lambda
    V_saved = V;
    lambda_saved = lambda;
    
    %% do voltage prediction to find predicted point (predicting voltage)
    [V_predicted, lambda_predicted, J] = cpf_predict(Ybus, ref, pv, pq, V, lambda, stepSize, 1, initQPratio, participation_i, flag_lambdaIncrease);
    
    %% do voltage correction to find corrected point
    [V, lambda, success, iterNum] = cpf_correctVoltage(baseMVA, bus, gen, Ybus, V_predicted, lambda_predicted, initQPratio, participation_i);

    %% calculate slope (dP/dLambda) at current point
	slope = min( abs( V-V_saved)./(lambda-lambda_saved));
	

    mean_error = mean( abs( V-V_predicted));
    error_order = log(mean_error/0.001);
    
    if abs(error_order) > 1 && mean(error)>0,
        newStepSize = stepSize * (1 + 0.8*(error_order < 0) - (0.8*error_order>0));
        newStepSize = max( min(newStepSize,maxStepSize),minStepSize); %clamp step size
        
		if verbose, fprintf('\t\tmean prediction error: %.15f. changed stepsize from %.5f to %.5f\n', mean_error, stepSize, newStepSize); end
        stepSize = newStepSize;
    end
   
    if lambda < 0 % lambda is less than 0, then stops CPF simulation
        if verbose > 0
            fprintf('\t[Info]:\tlambda is less than 0.\n\t\t\tCPF finished.\n');
        end
        break;
    end
    
    %% instead of using condition number as criteria for switching between
    %% modes...
    %%    if rcond(J) <= condNumThresh_Phase3 | success == false % Jacobian matrix is ill-conditioned, or correction step fails
    %% ...we use PV curve slopes as the criteria for switching modes
    
    if success
        if verbose > 2
            fprintf('\nVm_predicted\tVm_corrected\n');
            [[abs(V_predicted);lambda_predicted] [abs(V);lambda]]
        end

        %% record data
		logStepResults()
        pointCnt = pointCnt + 1;
        
    end
    
    
    if success == false % voltage correction step fails.
		V = V_saved;
        lambda = lambda_saved;
        if verbose > 0
			if ~isempty(strfind(lastwarn, 'singular'))
				fprintf('\t[Info]:\tMatrix is singular. Aborting Correction.\n');
				lastwarn('No error');
				break;
			else
				fprintf('\t[Info]:\tVoltage correction step fails..\n');
			end
        end
        break;        
    end
end


if success, %assuming we didn't fail out, try to solve for lambda = 0
    [V_predicted, lambda_predicted, ~] = cpf_predict(Ybus, ref, pv, pq, V, lambda_saved, lambda_saved, 1, initQPratio, participation_i, flag_lambdaIncrease);
    
    [V, lambda, success, ~] = cpf_correctVoltage(baseMVA, bus, gen, Ybus, V_predicted, lambda_predicted, initQPratio, participation_i);
    
    if success
        if verbose > 2
            fprintf('\nVm_predicted\tVm_corrected\n');
            [[abs(V_predicted);lambda_predicted] [abs(V);lambda]]
        end

        %% record data
        pointCnt = pointCnt + 1;
		logStepResults()
    end
end


max_lambda = max(lambda_corr);

pointCnt_Phase3 = pointCnt - pointCnt_Phase2 - pointCnt_Phase1; % collect number of points obtained at this phase



if verbose > 0, fprintf('\t[Info]:\t%d data points contained in phase 3.\n', pointCnt_Phase3); end
if shouldIPlotEverything, figure;	plot(lambda_corr, abs(V_corr)); end
et = etime(clock, t0);
if i == max_iter, fprintf('\t[Info] Max iterations hit.\n'); end


% if ~isempty(predicted_list) && ~isempty(corrected_list),
% 	%% combine the prediction and correction data in the sequence of appearance
% 	% NOTE: number of prediction data is one less than that of correction data
% 	predictedCnt = size(predicted_list, 2);
% 	combined_list(:, 1) = corrected_list(:, 1);
% 	for i = 1:predictedCnt
% 		combined_list(:, 2*i)     = predicted_list(:, i);
% 		combined_list(:, 2*i+1)   = corrected_list(:, i+1);
% 	end
% 
% 	%% convert back to original bus numbering & print results
% 	[bus, gen, branch] = int2ext(i2e, bus, gen, branch);
% 
% 	%% add bus number as the first column to the prediction, correction, and combined data list
% 	nb          = size(bus, 1);
% 	max_lambda  = max(corrected_list(nb+1, :));
% 	predicted_list = [[bus(:, BUS_I);0] predicted_list];
% 	corrected_list = [[bus(:, BUS_I);0] corrected_list];
% 	combined_list  = [[bus(:, BUS_I);0] combined_list];
% else
% 	combined_list = [];
% end


%% reorder according to exterior bus numbering
V_corr = V_corr(i2e(e2i),:);
Vpr = Vpr(i2e(e2i),:);

if nargout > 1, % LEGACY create predicted, corrected, combined lists
%     combined_list = zeros( size(V_corr,1) + 1, 1+size(V_corr,2) + size(Vpr,2)-1);
    predicted_list = [ [bus(:,1); 0] [Vpr; lambda_pr]];
    corrected_list = [ [bus(:,1); 0] [V_corr; lambda_corr]];
    
    combined_list = [bus(:,1); 0];
    combined_list(:,(1:size(V_corr,2))*2) = [V_corr; lambda_corr];
    combined_list(:,(1:size(Vpr,2)-1)*2+1) = [Vpr(:,2:end); lambda_pr(2:end)];
    
end

if verbose > 1
    Vm_corrected = abs(V_corr);
    Vm_predicted = abs(Vpr);
    Vm_corrected
    Vm_predicted
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
    
    function logStepResults()
        Vpr = [Vpr, V_predicted];
		lambda_pr = [lambda_pr, lambda_predicted];
		V_corr = [V_corr, V];
		lambda_corr = [lambda_corr, lambda];
		
		if shouldIPlotEverything, plotBusCurve(continuationBus); end
    end
    
    function plotBusCurve(bus)
        %use backwards order so that earlier points end up on top
        
        %plot phase 3
        plot(lambda_corr(1+i+j:1+i+j+k), abs(V_corr(bus, 1+i+j:1+i+j+k)), '.-b', 'markers',12); hold on;
        
        %plot phase 2
        plot(lambda_corr(1+i:1+i+j), abs(V_corr(bus, 1+i:1+i+j)), '.-g', 'markers', 12);
        
        %plot phase 1
        plot(lambda_corr(1:1+i), abs(V_corr(bus,1:i+1)), '.-b', 'markers', 12);
        
        
        %plot initial point
        plot(lambda_corr(1), abs(V_corr(bus,1)), '.-k', 'markers', 12);
             
        if 1+i+j+k < pointCnt,
           plot(lambda_corr(end-1:end), abs(V_corr(bus, end-1:end)),'.-k', 'markers', 12);
        end
%         fprintf('Points in 1, Phase 1, Phase 2, Phase 3: %d. Total Points: %d',1+i+j+k, pointCnt)
        
        scatter(lambda_pr, abs(Vpr(bus,:)),'r'); hold off;
        
    end


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
