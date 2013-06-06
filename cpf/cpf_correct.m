function [V, lambda, success, iter] = cpf_correct(caseData, state, V_predicted, lambda_predicted, continuationBus)


	switch state,
		case 1,
			[V, lambda, success, iter] = cpf_correctVoltage(caseData.baseMVA, caseData.bus, caseData.gen, caseData.Ybus, V_predicted, lambda_predicted, caseData.initQPratio, caseData.participation);
		case 2,
			 Vm_assigned = abs(V_predicted);
			[V, lambda, success, iter] = cpf_correctLambda(caseData.baseMVA, caseData.bus, caseData.gen, caseData.Ybus, Vm_assigned, V_predicted, lambda_predicted, caseData.initQPratio, caseData.participation, caseData.ref, caseData.pv, caseData.pq, continuationBus);
		case 3,
			[V, lambda, success, iter] = cpf_correctVoltage(caseData.baseMVA, caseData.bus, caseData.gen, caseData.Ybus, V_predicted, lambda_predicted, caseData.initQPratio, caseData.participation);
	end