function [V_predicted, lambda_predicted] = cpf_predict_voltage(V_corr, lambda_corr, lambda, sigma, ref,pv, pq,maxDegree)

    
    if nargin < 8,
        maxDegree = 6;
    end
    % degree of Lagrange polynomial is set by length of known data, too
		% high a degree causes instability

	
	%% set up indexing
	npv = length(pv);
	npq = length(pq);
	
	
	Vangles_pv = 1:npv;
	Vangles_pq = npv+1:npv+npq;
	Vmag_pq = npv+npq+1:npv+2*npq;
% 	lambdaIndex = npv+2*npq + 1;

	
	%% update lambda
	lambda_predicted = lambda + sigma;
	
	maxDegree = 10;
		
	 
	 
	V_corr =  V_corr(:, max(1,end-maxDegree+1):end);
	lambda_corr = lambda_corr(max(1,end-maxDegree+1):end);
		%shorten V_corr and lambda_corr to maxDegree points

		
		
		
	x_current = [ angle( V_corr([pv;pq], :) );  abs( V_corr(pq, : ))];
		
	x_predicted = lagrangepoly(lambda_predicted, lambda_corr, x_current);



	%% convert variable set to voltage form
	V_predicted(ref, 1) = V_corr(ref,end); %reference bus voltage passes through
	V_predicted(pv, 1) = abs(V_corr(pv, end)).* exp(sqrt(-1) * x_predicted(Vangles_pv) ); %apply new angle to 
	V_predicted(pq, 1) = x_predicted(Vmag_pq).* exp(sqrt(-1) * x_predicted(Vangles_pq) );
% 	lambda_predicted = x_predicted(lambdaIndex);

end



