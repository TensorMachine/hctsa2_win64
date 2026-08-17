function y = BF_idwt(a,d,wname,lx)
% BF_idwt   Single-level inverse DWT (toolbox-free idwt()).
% lx is the desired output length; the central lx samples are kept.
[~,~,Lo_R,Hi_R] = BF_wfilters(wname);
if nargin < 4 || isempty(lx)
    lx = 2*numel(a) - numel(Lo_R) + 2;
end
y = upsconv(a,Lo_R,lx) + upsconv(d,Hi_R,lx);
end

function z = upsconv(c,F,lx)
c = c(:);
z = zeros(2*numel(c)-1,1);   % dyadup: zeros interleaved
z(1:2:end) = c;
z = conv(z,F(:));
s = floor((numel(z)-lx)/2);  % wkeep, centred
z = z(s+1:s+lx);
end
