function [max_lambda, predicted_list, corrected_list, combined_list, success, et] = cpf(casedata, participation, sigmaForLambda, sigmaForVoltage, verbose, plotting)
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

GLOBAL_CONTBUS = 0; %for use by 'plotBusCurve'

lastwarn('No Warning');
%% assign default parameters
if nargin < 3
    sigmaForLambda = 0.1;       % stepsize for lambda
    sigmaForVoltage = 0.025;    % stepsize for voltage
end

if nargin < 5, verbose = 0; end
if nargin < 6, shouldIPlotEverything = false; else shouldIPlotEverything = plotting; end


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



if nargin < 2 %no participation factors so keep the current load profile
	participation = busE(:,PD)./sum(busE(:,PD));
else
	
	participation = participation(:);%(:) forces column vector

	if length(participation) ~= numBuses, %improper number of participations given
		if length(participation) == 1 && participation > 0,%assume bus number is specified instead
			participation = (1:numBuses)'==participation;
		else
			if verbose, fprintf('\t[Info]\tParticipation Factors improperly specified.\n\t\t\tKeeping Current Loading Profile.\n'); end
			participation = busE(:,PD)./sum(busE(:,PD));
		end
	end
end

%i2e is simply the bus ids stored in casedata.bus(:,1), so V_corr can be
%outputted as is with rows corresponding to entries in casedata.bus
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
    mError = addCause(mError, MException('MATPOWER:makeYBus', 'Ybus is singular'));
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





%%------------------------------------------------
% do cpf prediction-correction iterations
%%------------------------------------------------
t0 = clock;

nPoints = 0;
% V_pr=[];
% lambda_pr = [];
% V_corr = [];
% lambda_corr = [];

V_pr = zeros(size(bus,1), 400);
V_corr = zeros(size(bus,1),400);
lambda_pr = zeros(1,400);
lambda_corr = zeros(1,400);
stepSizes = zeros(1,400);





%% do voltage correction (ie, power flow) to get initial voltage profile
lambda_predicted = lambda;
V_predicted = V;
[V, lambda, success, ~] = cpf_correctVoltage(baseMVA, bus, gen, Ybus, V_predicted, lambda_predicted, initQPratio, participation_i);


if any(isnan(V))
    err = MException('CPF:correctVoltageError', 'Generating initial voltage profile');
    err = addCause(mError,MException('CPF:correctVoltageError', ['NaN bus voltage at ', mat2str(i2e(isnan(V)))]));    
    throw(err);
end

stepSize = 1;
logStepResults();
nPoints = nPoints + 1;









%% --- Start Phase 1: voltage prediction-correction (lambda increasing)
if verbose > 0, fprintf('Start Phase 1: voltage prediction-correction (lambda increasing).\n'); end

lagrange_order = 6;

%parametrize step size for Phase 1
minStepSize = 0.01;
maxStepSize = 1;

stepSize = sigmaForLambda;
stepSize = min(max(stepSize, minStepSize),maxStepSize);

finished = false;
continuationBus = pq(1);

function y= mean_log(x)
    y = log(x./mean(x));
end

i = 0; j=0; k=0; %initialize counters for each phase to zero
while i < max_iter    
    i = i + 1; % update iteration counter
    
    % save good data
    V_saved = V;
    lambda_saved = lambda;
    
	if nPoints<4 || any( abs(mean_log(stepSizes(end-lagrange_order:end)) > 1)), % do voltage prediction to find predicted point (predicting voltage)
		[V_predicted, lambda_predicted, ~] = cpf_predict(Ybus, ref, pv, pq, V, lambda, stepSize, 1, initQPratio, participation_i, flag_lambdaIncrease);
    else %if we have enough points, use lagrange polynomial
		[V_predicted, lambda_predicted] = cpf_predict_voltage(V_corr(:,1:nPoints), lambda_corr(1:nPoints), lambda, stepSize, ref, pv, pq, lagrange_order);
    end
    
           
    %% check prediction to make sure step is not too big so as to cause non-convergence    
    error_predicted = max(abs(V_predicted- V));
    if error_predicted > maxStepSize && ~success %-> this is inappropriate since 'success' would be coming from previous correction step
        newStepSize = 0.8*stepSize; %cut down the step size to reduce the prediction error
        if newStepSize > minStepSize,
            if verbose, fprintf('\t\tPrediction step too large (voltage change of %.4f). Step Size reduced from %.5f to %.5f\n', error_predicted, stepSize, newStepSize); end
            stepSize = newStepSize;
            i = i-1;            
            continue;
        end       
    end
    
    %% do voltage correction to find corrected point
    [V, lambda, success, ~] = cpf_correctVoltage(baseMVA, bus, gen, Ybus, V_predicted, lambda_predicted, initQPratio, participation_i);
	
    % if voltage correction fails, reduce step size and try again
    if success == false  && stepSize > minStepSize,
        newStepSize = stepSize * 0.3;
        if newStepSize > minStepSize,	
            if verbose, fprintf('\t\tCorrection step didnt converge; changed stepsize from %.5f to: %.5f\n', stepSize, newStepSize); end
            stepSize = newStepSize;
            i = i-1;
            V = V_saved;
            lambda = lambda_saved;
            continue;
        end
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
        logStepResults();        
		if shouldIPlotEverything, plotBusCurve(continuationBus); end
        nPoints = nPoints + 1;
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
            i = i-1;
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
if verbose > 0
    fprintf('\t[Info]:\t%d data points contained in phase 1.\n', i);
end








%% --- Switch to Phase 2: lambda prediction-correction (voltage decreasing)
if verbose > 0
    fprintf('Switch to Phase 2: lambda prediction-correction (voltage decreasing).\n');
end

busWhiteList = true(size(bus,1),1);

maxStepSize = 0.1;
minStepSize = 0.000001;
stepSize = sigmaForVoltage;
stepSize = min(max(stepSize, minStepSize),maxStepSize);
j = 0;
while j < max_iter && ~finished
    %% update iteration counter
    j = j + 1;

    % save good data
    V_saved = V;
    lambda_saved = lambda;
    
    
    %% do lambda prediction to find predicted point (predicting lambda)
	if nPoints<4 || any( abs(mean_log(stepSizes(end-lagrange_order:end)) > 1)),
        [V_predicted, lambda_predicted, J] = cpf_predict(Ybus, ref, pv, pq, V, lambda, stepSize, [2, continuationBus], initQPratio, participation_i,flag_lambdaIncrease);
    else
        [V_predicted, lambda_predicted] = cpf_predict_lambda(V_corr(:,1:nPoints), lambda_corr(1:nPoints), lambda, stepSize, continuationBus, ref, pv, pq, 4);
    end
   
    
    
    %% do lambda correction to find corrected point
    Vm_assigned = abs(V_predicted);
    
  
	
	[V, lambda, success, ~] = cpf_correctLambda(baseMVA, bus, gen, Ybus, Vm_assigned, V_predicted, lambda_predicted, initQPratio, participation_i, ref, pv, pq, continuationBus);
% 	
%     fprintf('Lambda error: %d
    
    if abs(lambda - lambda_saved) > maxStepSize, 
        %reject the new sample if lambda step is too big - this can happen
        %for certain buses. If it does happen, we can try reducing the step
        %size, discarding the sample and running the step again.
        
        if abs(lambda-lambda_saved) > 1 && busWhiteList(continuationBus) && stepSize < 0.0001, 
            % if lambda step is very large, there may be an issue with
            % using this bus as a continuation bus; therefore,
            % blacklist it and try a different bus
            busWhiteList(continuationBus) = false;
            if verbose, fprintf('\t\tBus %d blacklisted\n', continuationBus); end
            continuationBus = pickBus();
            
            V  = V_saved;
            lambda = lambda_saved;
            j = j-1;
            continue;
        end
        
        % ...otherwise, reduce the step size, discard and try again
        newStepSize = stepSize * 0.2;
        if newStepSize > minStepSize,
            if verbose, fprintf('continuationBus = %d', continuationBus); end
            if verbose, fprintf('\t\tLambda step too big; lambda step = %3f). Step size reduced from %.6f to %.6f\n', lambda-lambda_saved, stepSize, newStepSize); end
            
             
            stepSize = newStepSize;
            V = V_saved;
            lambda = lambda_saved;
            j = j-1;
            continue;
        end
         
    end
    
    
    
    
    
    %Here we check the change in Voltage if correction did not converge; if
    %the step is larger than the minimum then we can reduce the step size,
    %discard the sample and try again.
    mean_step = mean( abs(V_predicted-V_saved));    
	prediction_error = mean(abs(V-V_predicted)./abs(V));
    error_order = log(prediction_error/0.00001);
    
    if ( mean_step > 0.00001 && ~success) % if we jumped too far and correction step didn't converge
        newStepSize = stepSize * 0.4;
        if newStepSize > minStepSize, %if we are not below min step-size threshold go back and try again with new stepSize
            if verbose, fprintf('\t\tDid not converge; voltage step: %f pu. Step Size reduced from %.6f to %.6f\n', mean_step, stepSize, newStepSize); end
            stepSize = newStepSize;
            V = V_saved;
            lambda= lambda_saved;
            j =  j-1;
            continue;
        end
    end
    
    
    
    if abs(error_order) > 1.5 && prediction_error>0, 
        %if we havent just dropped the stepSize, consider changing it to
        %get a better error outcome. this allows us to increase our steps
        %to go faster or reduce our steps to avoid non-convergence, which
        %is time consuming.
        newStepSize = stepSize * (1 + 0.2*(error_order < 1) - 0.2*(error_order>1));
        newStepSize = max( min(newStepSize,maxStepSize),minStepSize); %clamp step size
        
		if verbose && newStepSize ~= stepSize, fprintf('\t\tAdjusting step size from %.6f to %.6f; mean prediction error: %.15f.\n', stepSize, newStepSize, prediction_error); end
        if newStepSize < 0.000001,
            keyboard
        end
        stepSize = newStepSize;
    end
	
    
	
    continuationBus = pickBus();

    
    if success %if correction step converged, log values and do verbosity
		logStepResults();        
		if shouldIPlotEverything, plotBusCurve(continuationBus); end
        nPoints = nPoints + 1;
    end
    
    
    
    % instead of using condition number as criteria for switching between
    % modes...
    %    if rcond(J) >= condNumThresh_Phase2 | success == false % Jacobian matrix is good-conditioned, or correction step fails
    % ...we use PV curve slopes as the criteria for switching modes:
    if ~success || (slope < 0 && slope > -slopeThresh_Phase2),
        if ~success,
            % restore good data
            V = V_saved;
            lambda = lambda_saved;
            j = j-1;
        end
        
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
if verbose > 0
    fprintf('\t[Info]:\t%d data points contained in phase 2.\n', j);
end


function continuationBus = pickBus()
    %how to pick the continuation bus during Phase 2:
    
    % 1. calculate slope (dP/dLambda) at current point
    mSlopes = abs(V-V_saved)./(lambda-lambda_saved);
    
    % 2. check if we have passed the peak of the PV curve
    if flag_lambdaIncrease && any(mSlopes < 0), flag_lambdaIncrease = false; end
    
    % 3. from pq buses that have not been blacklisted, pick the one with
    % the fastest changing variable
%     [~,continuationBusPQ] = max((-1)^(flag_lambdaIncrease+1) *mSlopes(pq(pqWhiteList))); %limit to only PQ busses
    
    mPQ = pq(busWhiteList(pq));
    if flag_lambdaIncrease,
        [~, ind] = max( mSlopes(mPQ));
    else
        [~, ind] = max( mSlopes(mPQ));
%         [~,ind] = mymedian(mSlopes(mPQ));
    end    
    continuationBus = mPQ(ind);
    
    slope = mSlopes(continuationBus);
end






%% --- Switch to Phase 3: voltage prediction-correction (lambda decreasing)
if verbose > 0
    fprintf('Switch to Phase 3: voltage prediction-correction (lambda decreasing).\n');
end
% set lambda to be decreasing
flag_lambdaIncrease = false; 


%set step size for Phase 3
minStepSize = 0.9*stepSize;
maxStepSize = 2;
stepSize = min(max(stepSize, minStepSize),maxStepSize);

k = 0;
while k < max_iter && ~finished
    %% update iteration counter
    k = k + 1;
    
    
    %% store V and lambda
    V_saved = V;
    lambda_saved = lambda;
    
    
    %% do voltage prediction to find predicted point (predicting voltage)
    [V_predicted, lambda_predicted, ~] = cpf_predict(Ybus, ref, pv, pq, V, lambda, stepSize, 1, initQPratio, participation_i, flag_lambdaIncrease);

    %% do voltage correction to find corrected point
    [V, lambda, success, ~] = cpf_correctVoltage(baseMVA, bus, gen, Ybus, V_predicted, lambda_predicted, initQPratio, participation_i);
	
    mean_step = mean( abs(V-V_saved));        
    if ( mean_step > 0.0001 && ~success) % if we jumped too far and correction step didn't converge
        newStepSize = stepSize * 0.4;
        if newStepSize > minStepSize, %if we are not below min step-size threshold go back and try again with new stepSize
            if verbose, fprintf('\t\tDid not converge; voltage step: %f pu. Step Size reduced from %.5f to %.5f\n', mean_step, stepSize, newStepSize); end
            stepSize = newStepSize;
            V = V_saved;
            lambda= lambda_saved;
            k =  k-1;
            continue;
        end
    end
    prediction_error = mean( abs( V-V_predicted));
    error_order = log(prediction_error/0.001);
    
    if abs(error_order) > 1 && prediction_error>0,
        newStepSize = stepSize - 0.03*log(prediction_error/0.001); %adjust step size

%         newStepSize = stepSize * (1 + 0.8*(error_order < 1) - 0.8*(error_order>1));
        newStepSize = max( min(newStepSize,maxStepSize),minStepSize); %clamp step size
        
   		if verbose && newStepSize ~= stepSize, fprintf('\t\tAdjusting step size from %.5f to %.5f; mean prediction error: %.15f.\n', stepSize, newStepSize, prediction_error); end
        
% 		if verbose, fprintf('\t\tmean prediction error: %.15f. changed stepsize from %.5f to %.5f\n', mean_error, stepSize, newStepSize); end
        stepSize = newStepSize;
    end
   
    if lambda < 0 % lambda is less than 0, then stop CPF simulation
        if verbose > 0, fprintf('\t[Info]:\tlambda is less than 0.\n\t\t\tCPF finished.\n'); end
        k = k-1;
        break;
    end
    
    
    if success,
        logStepResults()        
		if shouldIPlotEverything, plotBusCurve(continuationBus); end
        nPoints = nPoints + 1;        
    end
    
    
    if ~success % voltage correction step fails.
		V = V_saved;
        lambda = lambda_saved;
        k = k-1;
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
if verbose > 0, fprintf('\t[Info]:\t%d data points contained in phase 3.\n', k); end





%% Get the last point (Lambda == 0)
if success, %assuming we didn't fail out, try to solve for lambda = 0
    [V_predicted, lambda_predicted, ~] = cpf_predict(Ybus, ref, pv, pq, V, lambda_saved, lambda_saved, 1, initQPratio, participation_i, flag_lambdaIncrease);
    [V, lambda, success, ~] = cpf_correctVoltage(baseMVA, bus, gen, Ybus, V_predicted, lambda_predicted, initQPratio, participation_i);
    
    if success,        
		logStepResults()
        if shouldIPlotEverything, plotBusCurve(continuationBus); end
        
        nPoints = nPoints + 1;
    end
end

if verbose > 0, fprintf('\t[Info]:\t%d data points total.\n', nPoints); end

V_corr = V_corr(:,1:nPoints);
lambda_corr = lambda_corr(1:nPoints);
V_pr = V_pr(:,1:nPoints);
lambda_pr = lambda_pr(1:nPoints);

max_lambda = max(lambda_corr);




if shouldIPlotEverything, 
    plotBusCurve(continuationBus);
    figure;	
    
    hold on;      
        plot(lambda_corr, abs(V_corr(pq,:))); 
        maxL =plot([max_lambda, max_lambda], ylim,'LineStyle', '--','Color',[0.8,0.8,0.8]);  
%         mText = text(max_lambda*0.85, 0.1, sprintf('Lambda: %3.2f',max_lambda), 'Color', [0.7,0.7,0.7]);
        
        ticks = get(gca, 'XTick');
        ticks = ticks(abs(ticks - ml) > 0.5);
        ticks = sort(unique([ticks round(ml*1000)/1000]));
        set(gca, 'XTick', ticks);
    
        uistack(maxL, 'bottom');    
%         uistack(mText, 'bottom');
    hold off;
    
    
    title('Power-Voltage curves for all PQ buses.');
    ylabel('Voltage (p.u.)')
    xlabel('Power (lambda load scaling factor)');
    
%     ylims = ylim;
    
end
et = etime(clock, t0);


%% reorder according to exterior bus numbering

if nargout > 1, % LEGACY create predicted, corrected, combined lists
%     combined_list = zeros( size(V_corr,1) + 1, 1+size(V_corr,2) + size(Vpr,2)-1);
    predicted_list = [ [bus(:,1); 0] [V_pr; lambda_pr]];
    corrected_list = [ [bus(:,1); 0] [V_corr; lambda_corr]];
    
    combined_list = [bus(:,1); 0];
    combined_list(:,(1:size(V_corr,2))*2) = [V_corr; lambda_corr];
    combined_list(:,(1:size(V_pr,2)-1)*2+1) = [V_pr(:,2:end); lambda_pr(2:end)];
    
end


if nargout == 1,
	results.max_lambda = max_lambda;
	results.V_pr = V_pr;
	results.lambda_pr = lambda_pr;
	results.V_corr = V_corr;
	results.lambda_corr = lambda_corr;
	results.success = success;
	results.time = et;
	
	max_lambda = results; %return a single struct
end
    
    function logStepResults()
        
        V_pr(:,nPoints+1) = V_predicted;
        lambda_pr(nPoints+1) = lambda_predicted;
        V_corr(:,nPoints+1) = V;
        lambda_corr(:,nPoints+1) = lambda;
        stepSizes(:,nPoints+1) = stepSize;
    end
    
    function plotBusCurve(bus)
        
        if bus == GLOBAL_CONTBUS, %if its the same bus as last time, check resizing of window.
            xlims = xlim;
            ylims = ylim;

            xlims(2) = max(xlims(1) + (lambda_corr(nPoints)- xlims(1)) * 1.2, xlims(2));
            ylims(1) = min(ylims(2) - (ylims(2) - abs(V_corr(bus,nPoints)))*1.2, ylims(1));
        end
        %plot phase 3
        p3=plot(lambda_corr(1+i+j:1+i+j+k), abs(V_corr(bus, 1+i+j:1+i+j+k)), '.-b', 'markers',12); hold on;
        
        %plot phase 2
        p2=plot(lambda_corr(1+i:1+i+j), abs(V_corr(bus, 1+i:1+i+j)), '.-g', 'markers', 12);
        
        %plot phase 1
        p1=plot(lambda_corr(1:1+i), abs(V_corr(bus,1:i+1)), '.-b', 'markers', 12);
        
        
        %plot initial point
        st=plot(lambda_corr(1), abs(V_corr(bus,1)), '.-k', 'markers', 12);
             
        if 1+i+j+k < nPoints,
            en=plot(lambda_corr(end-1:end), abs(V_corr(bus, end-1:end)),'.-k', 'markers', 12);
            
            scatter(lambda_pr(nPoints), abs(V_pr(bus, nPoints)), 'r');
        end
%         fprintf('Points in 1, Phase 1, Phase 2, Phase 3: %d. Total Points: %d',1+i+j+k, pointCnt)
        
        pred=scatter(lambda_pr(1:1+i+j+k), abs(V_pr(bus,1:1+i+j+k)),'r'); hold off;
        
        
        title(sprintf('Power-Voltage curve for bus %d.', continuationBus));
        ylabel('Voltage (p.u.)')
        xlabel('Power (lambda load scaling factor)');
        ml = max(lambda_corr);
        
        ticks = get(gca, 'XTick');
        ticks = ticks(abs(ticks - ml) > 0.25);
        ticks = sort(unique([ticks round(ml*1000)/1000]));
        set(gca, 'XTick', ticks);
        hold on;
            mLine=plot( [ml, ml], ylim, 'LineStyle', '--', 'Color', [0.8,0.8,0.8]);
            uistack(mLine, 'bottom');
        hold off;
        legend([pred, p1,p2, mLine], {'Predicted Values', 'Lambda Continuation', 'Voltage Continuation', 'Max Lambda'})
        
        if bus == GLOBAL_CONTBUS,
            xlim(xlims);
            ylim(ylims);
            
        end
        
        GLOBAL_CONTBUS = bus;
        a = 1;
    end

function [med, idx] =mymedian(x)
% mymedian    Calculate the median value of an array and find its index.
%
% This function can be used to calculate the median value of an array.
% Unlike the built in median function, it returns the index where the
% median value occurs. In cases where the array does not contain its
% median, such as [1,2,3,4] or [1,3,2,4], the index of the first occuring
% adjacent point will be returned; in both of the above examples the median
% will be 2.5 and the index will be 2.

    assert(isvector(x));
    
    
    med = median(x);
    
    [~, idx] = min( abs( x - med));
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
