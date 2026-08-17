function [Lo_D,Hi_D,Lo_R,Hi_R] = BF_wfilters(wname)
% BF_wfilters   Toolbox-free replacement for the Wavelet Toolbox's wfilters().
%
% Only the two wavelets hctsa actually uses are provided: db3 and sym2.
% Coefficients are the standard orthonormal scaling filters; the quadrature
% mirror and reversal conventions match MATLAB's orthfilt().

switch lower(wname)
case 'db3'
    Lo_R = [ 0.3326705529500825,  0.8068915093110924,  0.4598775021184914, ...
            -0.1350110200102546, -0.0854412738820267,  0.0352262918857095];
case {'haar','db1'}
    % Haar / db1. Used by wfbmesti's third estimate.
    Lo_R = [1 1]/sqrt(2);
case 'sym5'
    % Symlet-5 reconstruction low-pass filter. Needed by wfbmesti, which takes
    % the sym5 HIGH-PASS decomposition filter (the second wfilters output).
    Lo_R = [ 0.0195388827352869, -0.0211018340249425, -0.1753280899084594, ...
             0.0166021057644804,  0.6339789634582119,  0.7234076904024206, ...
             0.1993975339773936, -0.0391342493023831,  0.0295194909257746, ...
             0.0273330683450780];
case {'sym2','db2'}
    Lo_R = [ 0.4829629131445341,  0.8365163037378079,  0.2241438680420134, ...
            -0.1294095225512604];
otherwise
    error('BF_wfilters:unsupported', ...
          ['Wavelet ''%s'' is not supported. hctsa 2.0 provides haar, db3, sym2 ' ...
           'and sym5, which are the wavelets used by the default feature ' ...
           'library.'],wname);
end

% Quadrature mirror: Hi_R(k) = (-1)^(k-1) * Lo_R(end-k+1)
Hi_R = Lo_R(end:-1:1) .* (-1).^(0:numel(Lo_R)-1);
Lo_D = Lo_R(end:-1:1);
Hi_D = Hi_R(end:-1:1);
end
