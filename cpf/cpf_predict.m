function [V_predicted, lambda_predicted, J] = cpf_predict(Ybus, ref, pv, pq, V, lambda, sigma, type_predict, initQPratio, loadvarloc, flag_lambdaIncrease)
%CPF_PREDICT  Do prediction in cpf.
%   [INPUT PARAMETERS]
%   type_predict: 1-predict voltage; 2-predict lambda
%   loadvarloc: (in internal bus numbering)
%   [OUTPUT PARAMETERS]
%   J: jacobian matrix for the given voltage profile (before prediction)
%   created by Rui Bo on 2007/11/12

%   MATPOWER
%   $Id: cpf_predict.m,v 1.4 2010/04/26 19:45:26 ray Exp $
%   by Rui Bo
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

participation = loadvarloc;

%% set up indexing
npv = length(pv);
npq = length(pq);

pv_bus = ~isempty(find(pv == loadvarloc));

%% form current variable set from given voltage
x_current = [ angle(V([pv;pq]));
              abs(V(pq));
              lambda];

%% evaluate Jacobian
[dSbus_dVm, dSbus_dVa] = dSbus_dV(Ybus, V);

j11 = real(dSbus_dVa([pv; pq], [pv; pq]));
j12 = real(dSbus_dVm([pv; pq], pq));
j21 = imag(dSbus_dVa(pq, [pv; pq]));
j22 = imag(dSbus_dVm(pq, pq));

J = [   j11 j12;
        j21 j22;    ];


% form K based on participation factors.
participation = participation ./sum(participation); %normalize
K = zeros(npv+2*npq,1);
% 
% K(pv) = participation(pv);
% K(npv+pq) = participation(pq);
% K(npv+npq+pq) = participation(pq).*initQPratio;
K(1:length(pv)) = participation(pv);
K(npv + (1:length(pq))) = -participation(pq);
K(npv+npq+(1:length(pq))) = -participation(pq) .* initQPratio(pq);

% %% form K
% K = zeros(npv+2*npq, 1);
% if pv_bus % pv bus
%     K(find(pv == loadvarloc)) = -1;                         % corresponding to deltaP
% else % pq bus
%     K(npv + find(pq == loadvarloc)) = -1;                   % corresponding to deltaP
%     K(npv + npq + find(pq == loadvarloc)) = -initQPratio;   % corresponding to deltaQ
% end

%% form e

e = zeros(1, npv+2*npq+1);
if type_predict(1) == 1 % predict voltage
    if flag_lambdaIncrease == true
        e(npv+2*npq+1) = 1; % dLambda = 1
    else
        e(npv+2*npq+1) = -1; % dLambda = -1
    end
elseif type_predict(1) == 2 % predict lambda
    % [Anton] we have to discriminate between pv and pq bus because in the
	% original CPF, the changing load was used as the voltage continuation
	% parameter, so using pv bus would be nonsensical. now we use whichever
	% bus changes the most, so pv bus must be considered
	
	
	%each bus has an angle, plus all PQ busses have a Voltage
	continuationBus = type_predict(2);%% [Anton] used type_predict to pass in bus value I want for voltage continuation
	
	if any(pq == continuationBus),
		e(npv + npq + find(pq == continuationBus)) = -1;
	elseif any(pv==continuationBus),
		e(find(pv==continuationBus)) = -1;
	end
    %e(npv+npq+find(pq == loadvarloc)) = -1; % dVm = -1
else
    fprintf('Error: unknow ''type_predict''.\n');
    pause
end

% form of e is expected to be [ delta * (#pv buses + #pq buses), v* #pq
% buses) + 1]

%% form b
b = zeros(npv+2*npq+1, 1);
b(npv+2*npq+1) = 1;

%% form augmented Jacobian
%NOTE: the use of '-J' instead of 'J' is due to that the definition of
%dP(,dQ) in the textbook is the negative of the definition in MATPOWER. In
%the textbook, dP=Pinj-Pbus; In MATPOWER, dP=Pbus-Pinj. Therefore, the
%Jacobians generated by the two definitions differ only in the sign.
augJ = [-J K;   
        e   ];

%% calculate predicted variable set
x_predicted = x_current + sigma*(augJ\b);

%% convert variable set to voltage form
V_predicted([ref], 1) = V([ref]);
V_predicted([pv], 1) = abs(V([pv])).* exp(sqrt(-1) * x_predicted([1:npv]) );
V_predicted([pq], 1) = x_predicted([npv+npq+1:npv+2*npq]).* exp(sqrt(-1) * x_predicted([npv+1:npv+npq]) );
lambda_predicted = x_predicted(npv+2*npq+1);
