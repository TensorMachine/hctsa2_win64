function [C,L] = BF_wavedec(x,level,wname)
% BF_wavedec   Multilevel 1-D wavelet decomposition (toolbox-free wavedec()).
%
% C = [cA_level, cD_level, cD_level-1, ..., cD_1]
% L = [len(cA_level), len(cD_level), ..., len(cD_1), len(x)]

x = x(:);
C = [];
L = numel(x);
a = x;
for k = 1:level
    [a,d] = BF_dwt(a,wname);
    C = [d(:); C];
    L = [numel(d); L];
end
C = [a(:); C];
L = [numel(a); L];
C = C(:)';
L = L(:)';
end
