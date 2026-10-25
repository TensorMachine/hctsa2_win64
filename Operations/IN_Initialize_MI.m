function miCalc = IN_Initialize_MI(estMethod, extraParam, addNoise, y)
% IN_Initialize_MI  Initialize a mutual information calculator (JIDT's interface).
%
% HCTSA2 2.4.0: returns a BF_MutualInfoCalculator; the JIDT jar is no longer shipped.
%
% ---INPUTS:
%
% estMethod: the estimation method used to compute the mutual information:
%           (*) 'gaussian'
%           (*) 'kernel' (JIDT's; not available since hctsa2 2.4.0)
%           (*) 'kraskov1'
%           (*) 'kraskov2'
%
% y: (optional) the data about to be analyzed, used only to check for a high
%    proportion of repeated values when addNoise is false (see below).
%
% cf. Kraskov, A., Stoegbauer, H., Grassberger, P., Estimating mutual
% information: http://dx.doi.org/10.1103/PhysRevE.69.066138

% ------------------------------------------------------------------------------
% Copyright (C) 2013-2026, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite the following two papers:
%
% (1) B.D. Fulcher and N.S. Jones, "hctsa: A Computational Framework for Automated
% Time-Series Phenotyping Using Massive Feature Extraction, Cell Systems 5: 527 (2017).
% DOI: 10.1016/j.cels.2017.10.001
%
% (2) B.D. Fulcher, M.A. Little, N.S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013).
% DOI: 10.1098/rsif.2013.0048
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program. If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

if nargin < 2
	extraParam = [];
end
if nargin < 3
	addNoise = false; % deterministic by default
end
if nargin < 4
	y = [];
end

% ------------------------------------------------------------------------------
% HCTSA2 2.4.0: the calculator is a BF_MutualInfoCalculator, a MATLAB class with
% JIDT's method names whose Kraskov estimator is the KSG_mi MEX (built from source
% by HCTSA2_BuildBinaries), instead of a class from the JIDT jar, which hctsa2 no
% longer ships. JIDT's kernel estimator has no replacement: no registered feature
% uses it.
switch estMethod
	case {'gaussian', 'kraskov1', 'kraskov2'}
	case 'kernel'
		error('IN_Initialize_MI:kernel', ['The ''kernel'' estimator came from JIDT, ' ...
			'which hctsa2 no longer ships (2.4.0); use ''kraskov1'', ''kraskov2'' or ''gaussian''.']);
	otherwise
		error('Unknown mutual information estimation method ''%s''', estMethod);
end
% ------------------------------------------------------------------------------

% Add neighest neighbor option for KSG estimator
k = 4; % use 4 nearest neighbors for KSG estimator as default
% (empirically the best of k=2/3/4 on ground-truth AR(1) and periodic
% test signals -- lowest RMSE/variance, no bias cost -- and matches
% what every registered hctsa feature already hard-codes)
if ismember(estMethod, {'kraskov1', 'kraskov2'}) && ~isempty(extraParam)
	if isnumeric(extraParam)
		warning(['Should set number of nearest neighbors to a string... ' ...
				 'I''ll do it for you this time though...'])
		extraParam = num2str(extraParam);
	end
	k = str2double(extraParam); % 4th input specifies number of nearest neighbors for KSG estimator
end

% Deterministic: no noise is added inside the estimator. Tie-breaking jitter for
% data with many repeated values is added by the caller, reproducibly: see
% BF_TieBreakNoise. (JIDT's own NOISE_LEVEL_TO_ADD drew from an unseeded Java RNG.)
if addNoise && ismember(estMethod, {'kraskov1', 'kraskov2'})
	error('IN_Initialize_MI:noise', ['Noise added inside the estimator is not ' ...
		'available (hctsa2 2.4.0); see BF_TieBreakNoise.']);
end

miCalc = BF_MutualInfoCalculator(estMethod, k);

% Specify a univariate calculation:
miCalc.initialise(1, 1);

end
