function [a,d] = BF_dwt(x,wname)
% BF_dwt   Single-level discrete wavelet transform (toolbox-free dwt()).
%
% Uses MATLAB's default half-point symmetric boundary extension ('sym').

[Lo_D,Hi_D] = BF_wfilters(wname);
x = x(:);
lf = numel(Lo_D);
n = numel(x);
if n < 1
    error('BF_dwt:emptyInput','The input signal is empty.');
end

% Half-point symmetric extension by lf-1 on each side, by index folding: for
% n < lf-1 the reflection repeats, as in wextend('sym') (HCTSA2 2.2.1).
lE = lf - 1;
m = mod((1-lE:n+lE) - 1,2*n);        % 0-based position within one period of 2n
ix = m + 1;
ix(m >= n) = 2*n - m(m >= n);        % reflected half of the period
xe = x(ix(:));

a = downsample2(conv(xe,Lo_D(:),'valid'));
d = downsample2(conv(xe,Hi_D(:),'valid'));
end

function y = downsample2(z)
y = z(2:2:end);
end
