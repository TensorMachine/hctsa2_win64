function [a,d] = BF_dwt(x,wname)
% BF_dwt   Single-level discrete wavelet transform (toolbox-free dwt()).
%
% Uses MATLAB's default half-point symmetric boundary extension ('sym').

[Lo_D,Hi_D] = BF_wfilters(wname);
x = x(:);
lf = numel(Lo_D);
n = numel(x);

% Half-point symmetric extension by lf-1 on each side
lE = lf - 1;
left  = x(min(lE,n):-1:1);
right = x(n:-1:max(1,n-lE+1));
if numel(left) < lE,  left  = [repmat(x(1),lE-numel(left),1);  left];  end
if numel(right) < lE, right = [right; repmat(x(end),lE-numel(right),1)]; end
xe = [left; x; right];

a = downsample2(conv(xe,Lo_D(:),'valid'));
d = downsample2(conv(xe,Hi_D(:),'valid'));
end

function y = downsample2(z)
y = z(2:2:end);
end
