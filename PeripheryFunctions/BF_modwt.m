function w = BF_modwt(x,wname,level)
% BF_modwt   Maximal overlap discrete wavelet transform (toolbox-free modwt()).
%
%   w = BF_modwt(x,wname,level)
%
% Returns a (level+1)-by-N matrix: rows 1..level are the wavelet coefficients
% at levels 1..level, the last row the scaling coefficients at the final level.
% Periodic boundary (MATLAB's default), which is the only one hctsa uses.
%
% The Percival & Walden pyramid algorithm (Wavelet Methods for Time Series
% Analysis, 2000, section 5.5) with the MODWT filters h/sqrt(2), g/sqrt(2),
% where h and g are the orthogonal wavelet's reconstruction high- and low-pass
% filters, and circular filtering of the previous level's scaling coefficients
% by the filters upsampled by 2^(j-1):
%
%   W(j,t) = sum_l h(l)/sqrt(2) * V(j-1, t - 2^(j-1)*l mod N)
%   V(j,t) = sum_l g(l)/sqrt(2) * V(j-1, t - 2^(j-1)*l mod N),   V(0,:) = x
%
% Matches MATLAB's modwt(x,'db3',level) and modwt(x,'sym4',level) to 1e-11
% (the size of the difference between MATLAB's stored db3 coefficients and
% the exact ones) on the eight series of the licensed run of 2026-09-23.

x = x(:)';
N = numel(x);
if nargin < 3 || isempty(level)
    level = floor(log2(N));
end
if level < 1 || level > floor(log2(N))
    error('BF_modwt:level','Level must be between 1 and floor(log2(N)) = %u.',floor(log2(N)));
end

[~,~,Lo_R,Hi_R] = BF_wfilters(wname);
h = Hi_R(:)'/sqrt(2);
g = Lo_R(:)'/sqrt(2);

w = zeros(level+1,N);
V = x;
for j = 1:level
    step = 2^(j-1);
    Wj = zeros(1,N);
    Vj = zeros(1,N);
    for l = 1:numel(h)
        Vs = circshift(V,[0, mod(step*(l-1),N)]);    % Vs(t) = V(t - step*(l-1))
        Wj = Wj + h(l)*Vs;
        Vj = Vj + g(l)*Vs;
    end
    w(j,:) = Wj;
    V = Vj;
end
w(level+1,:) = V;
end
