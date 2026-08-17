function [pacf,lags,bounds] = BF_parcorr(y,varargin)
% BF_parcorr   Sample partial autocorrelation (toolbox-free parcorr()).
%
% Method 'ols' (MATLAB's default): for each lag k, regress y on a constant and
% its first k lags, and take the coefficient on the kth lag.
% Method 'yule-walker': Levinson-Durbin recursion on the sample ACF.
%
% pacf(1) is lag 0 and equals 1.

y = y(:); y = y(isfinite(y));
N = numel(y);
numLags = min(20,N-1);
method = 'ols';
for k = 1:2:numel(varargin)-1
    switch lower(varargin{k})
    case 'numlags', numLags = varargin{k+1};
    case 'method',  method  = varargin{k+1};
    end
end

pacf = zeros(numLags+1,1);
pacf(1) = 1;

if strncmpi(method,'y',1)
    acf = BF_autocorr(y,'NumLags',numLags);
    phi = zeros(numLags,numLags);
    phi(1,1) = acf(2);
    pacf(2) = phi(1,1);
    for k = 2:numLags
        num = acf(k+1) - sum(phi(k-1,1:k-1)'.*acf(k:-1:2));
        den = 1 - sum(phi(k-1,1:k-1)'.*acf(2:k));
        phi(k,k) = num/den;
        phi(k,1:k-1) = phi(k-1,1:k-1) - phi(k,k)*phi(k-1,k-1:-1:1);
        pacf(k+1) = phi(k,k);
    end
else
    for k = 1:numLags
        yv = y(k+1:end);
        Xv = ones(numel(yv),1);
        for j = 1:k
            Xv = [Xv, y(k+1-j:end-j)];
        end
        b = Xv\yv;
        pacf(k+1) = b(end);
    end
end
lags = (0:numLags)';
bounds = [2;-2]/sqrt(N);
end
