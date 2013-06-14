function [V, lambda, converged, iterNum] = cpf_correctLambda(baseMVA, bus, gen, Ybus, Vm_assigned, V_predicted, lambda_predicted, initQPratio, participation, ref, pv, pq, continuationBus)
%CPF_CORRECTLAMBDA  Correct lambda in correction step near load point.
%   function: correct lambda(ie, real power of load) in cpf correction step
%   near the nose point. Use NR's method to solve the nonlinear equations
%   [INPUT PARAMETERS]
%   loadvarloc: (in internal bus numbering)
%   created by Rui Bo on 2007/11/12

%   MATPOWER
%   $Id: cpf_correctLambda.m,v 1.4 2010/04/26 19:45:26 ray Exp $
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
% 	[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
% 		VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
	[~, ~, ~, ~, ~, ~, PD, QD, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = idx_bus;


	%% options
	tolerance     = 1e-5; % mpopt(2);
	max_iters  = 100;  % mpopt(3);
	verbose = 0;    % mpopt(31);
	
	%% initialize
	V = V_predicted;
	lambda = lambda_predicted;
	Va = angle(V);
	Vm = abs(V);

	%% set up indexing for updating V
	npv = length(pv);
	npq = length(pq);
	
	Vangles_pv = 1:npv;
	Vangles_pq = npv+1:npv+npq;
	Vmag_pq = npv+npq+1:npv+2*npq;
	lambdaIndex = npv+2*npq + 1;
		
	%% update bus and calculate power flow value at current V
	[bus, F] = updatePF(bus);
	
	%% do Newton iterations
	i = 0;
	converged = false;
	while (~converged && i < max_iters)
		i = i + 1;

		%% evaluate Jacobian
		J = getJ(Ybus, V, pv, pq);

		%% form augmented Jacobian with V as continuation parameter
		delF_dLambda = [-participation(pv); -participation(pq); -participation(pq) .* initQPratio(pq)];
 		delVm = [zeros(1, npv+npq) -participation(pq)' 0];
			%dV/dangle = 0,  dV/dV = participations, dV/dlambda = 0;
		
		delVm = zeros(1, npv+npq*2 + 1);
% 		delVm(npv+npq+find(pq == continuationBus)) = -1;
		delVm(npv+find(pq == continuationBus)) = -1;
		delVm(pv == continuationBus) = -1;
		
		augJ = [ J delF_dLambda;
				  delVm          ];

		%% compute update step
		dx = -(augJ \ F);

		if ~isempty(strfind(lastwarn, 'singular')), 
			lastwarn('No error');
			break;
			converged = false;
		end

		%% update voltage. 
		Va( [pv;pq] ) = Va( [pv;pq] ) + dx( [Vangles_pv Vangles_pq]);
		Vm(pq) = Vm(pq) + dx(Vmag_pq);
		lambda = lambda + dx(lambdaIndex);
			% NOTE: voltage magnitude of pv buses, voltage magnitude
			% and angle of reference bus are not updated, so they keep as constants
			% (ie, the value as in the initial guess)

		V = Vm .* exp(1i * Va); % NOTE: angle is in radians in pf solver, but in degree in case data
		Vm = abs(V);            %% update Vm and Va again in case
		Va = angle(V);          %% we wrapped around with a negative Vm
 		
		
		[bus, F] = updatePF(bus);
		
		%% check for convergence
		normF = norm(F, inf);
		if verbose > 1,	fprintf('\niteration [%3d]\t\tnorm of mismatch: %10.3e', i, normF);	end
		converged = normF < tolerance;

	end

	iterNum = i;
	
% 	error(lastwarn);
	
	if ~isempty(strfind(lastwarn, 'singular')), 
		lastwarn('No error');
		converged = false;
	end

	
	
	function [bus, F] = updatePF(bus)
		%% set load as lambda indicates
		bus(:,PD) = bus(:,PD) .*(participation == 0) + participation.*lambda.*baseMVA;
		bus(:,QD) = bus(:,PD) .* initQPratio;

		%% compute complex bus power injections (generation - load)
		SbusInj = makeSbus(baseMVA, bus, gen);

		%% evalute F(x0)
		F = Feval(V,Vm_assigned, SbusInj, Ybus, pv, pq, continuationBus);
	end
end




function F = Feval(V, Vm_assigned, SbusInj, Ybus, pv, pq, continuationBus)
% F = Feval(V, Vm_assigned, SbusInj, Ybus, pv, pq, continuationBus)
%
% Feval   evaluate power flow equations
%
% Evaluates the power flow equations  S = V x V* / X
	
	%% calculate mismatch
	mis = SbusInj -  (  V .* conj(Ybus * V)  );
	
	%% order mismatches according to formulation
	F = [   real(mis(pv));
			real(mis(pq));
			imag(mis(pq));
			abs(V(continuationBus)) - Vm_assigned(continuationBus);   ];
end

function J = getJ(Ybus, V, pv, pq)
% J = getJ(Ybus, V, pv, pq)
%
% getJ   get jacobian of power flow equations
%
%  This function gets the Jaciobian of the power flow equations at 
	[dSbus_dVm, dSbus_dVa] = dSbus_dV(Ybus, V);

	j11 = real(dSbus_dVa([pv; pq], [pv; pq]));
	j12 = real(dSbus_dVm([pv; pq], pq));
	j21 = imag(dSbus_dVa(pq, [pv; pq]));
	j22 = imag(dSbus_dVm(pq, pq));

	J = -[   j11 j12;
			j21 j22;    ];
	%% form augmented Jacobian
	%NOTE: the use of '-J' instead of 'J' is due to that the definition of
	%dP(,dQ) in the textbook is the negative of the definition in MATPOWER. In
	%the textbook, dP=Pinj-Pbus; In MATPOWER, dP=Pbus-Pinj. Therefore, the
	%Jacobians generated by the two definitions differ only in the sign.
end
