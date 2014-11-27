function   Lx = lagrangepoly(x, xs, ys)
	% lagrangepoly	Extrapolate/Interpolate functions using lagrange polynomial.
	%
    % Lx = lagrangepoly(x,xs,ys);
    %
	% This function uses the lagrange polynomial formulation to perform interpolation
	% or extrapolation via polynomial. It works by returning y = L(x), where L represents the lagrange 
	% polynomial derived from x,xs and ys. The order of the polynomial is
	% determined by the number of samples in the input.
	%
	% inputs:
	%  x: x-axis value for which you are seeking y value inter/extra-polation.
	%  xs: x values corresponding to input functional values, in columns
	%  ys: y data from which to extrapolate the output y=L(x), where L is the Lagrange Polynomial.
	%      each column of ys is considered as a separate function, with one sample per column; thus xs and ys should 
	%      have the same  number of columns, and y will have the same numer of rows as ys.

	nSamples = length(xs);

	xxm = x-xs;


	js = 1:nSamples;

	Lj = zeros(nSamples,1);
	for j = 1:nSamples,
		Lj(j) = prod(xxm(:,js~=j),2) ./  prod( xs(j)- xs(:, js~=j));


	end

	Lx = ys * Lj;
    
%     figure; plot(xs, ys(21,:), 'b.-');
%     hold on; plot(x, Lx(21), 'ro'); hold off;
    

end
% 
% function  y = lagrangepoly(x, xs, ys)
% 
% 	ls  = zeros(1, length(xs));
% 	
% 	for j = 1:length(xs),
% 		ls(j) = ys(j) * lagrange_coefficient(x,j,xs,ys);
% 	end
% 	
% 	y = sum(ls);
% 		
% 		
% 
% 
% 
% end
% 
% function lj =  lagrange_coefficient(x,j, xs, ys)
% 	
% 	notj = 1:length(xs); notj = notj(notj ~= j);
% 	
% 	top = x - xs(notj);
% 	bottom = xs(j) - xs(notj);
% 	
% 	lj = prod(top)/prod(bottom);
% end