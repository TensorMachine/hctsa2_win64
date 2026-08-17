function Hest = BF_wfbmesti(x)
% BF_wfbmesti   Fractal index estimates (wfbmesti()).
%
% Toolbox-free replacement, written against MathWorks' wfbmesti.m source. Three
% estimates of the fractal index H of a signal taken to come from a fractional
% Brownian motion:
%
%   1. second-order discrete derivative, filters [1 -2 1] and [1 0 -2 0 1]
%   2. the same construction but with the sym5 **high-pass decomposition**
%      filter, against that filter upsampled by 2
%   3. linear regression of log2(detail variance) on level, using **haar** and
%      the robust (MAD-based) deviation estimate from wnoisest
%
% Estimates 1 and 2 operate on cumsum(diff(x)); estimate 3 uses x directly.

x = x(:);
x = x(isfinite(x));

% Estimates 1 and 2 work on the integrated first difference
y = diff(x);
y = cumsum(y(:)');
n = length(y);

Hest = zeros(1,3);

%-------------------------------------------------------------------------------
%% 1. second-order discrete derivative
%-------------------------------------------------------------------------------
b1 = [1 -2 1];
b2 = [1  0 -2 0 1];
y1 = filter(b1,1,y);  y1 = y1(length(b1):n);
y2 = filter(b2,1,y);  y2 = y2(length(b2):n);
Hest(1) = 0.5*log2(mean(y2.^2)/mean(y1.^2));

%-------------------------------------------------------------------------------
%% 2. same, using the sym5 high-pass decomposition filter
%-------------------------------------------------------------------------------
% wfilters returns [Lo_D, Hi_D, Lo_R, Hi_R]; wfbmesti takes the SECOND output.
[~,c1] = BF_wfilters('sym5');
c1 = c1(:)';
c2 = [c1; zeros(1,length(c1))];
c2 = c2(:)';                       % c1 upsampled by 2
cy1 = filter(c1,1,y);  cy1 = cy1(length(c1):n);
cy2 = filter(c2,1,y);  cy2 = cy2(length(c2):n);
Hest(2) = 0.5*log2(mean(cy2.^2)/mean(cy1.^2));

%-------------------------------------------------------------------------------
%% 3. variance versus level (haar)
%-------------------------------------------------------------------------------
levdec = min(BF_wmaxlev(numel(x),'haar'),6);
if levdec < 1
    Hest(3) = NaN;
    return
end
[c,l] = BF_wavedec(x,levdec,'haar');
lvls = 1:levdec;
stdc = BF_wnoisest(c,l,lvls);
p = polyfit(lvls,log2(stdc.^2),1);
Hest(3) = (p(1)-1)/2;
end
