function drawPV(CPFresults, bussesToDraw)
%drawPV  Draw PV curves for specified buses.
%   [INPUT PARAMETERS]
%   corrected_list, combined_list: data points obtained from CPF solver
%   loadvarloc: load variation location(in external bus numbering). Single bus supported so far.
%   flag_combinedCurve: flag indicating if the prediction-correction curve will be drawn
%   busesToDraw: bus indices whose PV curve will be be drawn
%   created by Rui Bo on 2008/01/13

%   MATPOWER
%   $Id: drawPVcurves.m,v 1.6 2010/04/26 19:45:26 ray Exp $
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

%% assign default parameters

if nargin < 2,
	bussesToDraw = 1:size(CPFresults.V_corr, 1)';
end

bussesToDraw = bussesToDraw(bussesToDraw <= size(CPFresults.V_corr,1));


plot(CPFresults.lambda_corr, abs(CPFresults.V_corr(bussesToDraw, :)));

title('CPF curves');
xlabel('lambda'); ylabel('Bus voltage (p.u.)');

