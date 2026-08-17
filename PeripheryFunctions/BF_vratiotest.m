function [h,pValue,stat,cValue,ratio] = BF_vratiotest(y,varargin)
% BF_vratiotest   Lo-MacKinlay variance ratio test (toolbox-free vratiotest()).
%
% Y is treated as a log price series. Under the random walk null the variance of
% q-period differences is q times the variance of 1-period differences.
%
% Name/value pairs: 'period' (default 2), 'IID' (default false), 'alpha' (0.05).
%   IID true  -> homoskedastic null
%   IID false -> heteroskedasticity-robust null

y = y(:);
y = y(isfinite(y));
T = numel(y);

period = 2; iid = false; alpha = 0.05;
for k = 1:2:numel(varargin)-1
    switch lower(varargin{k})
    case 'period', period = varargin{k+1};
    case 'iid',    iid    = varargin{k+1};
    case 'alpha',  alpha  = varargin{k+1};
    end
end
period = period(:)'; iid = logical(iid(:))'; alpha = alpha(:)';
n = max([numel(period) numel(iid) numel(alpha)]);
if numel(period)==1, period = repmat(period,1,n); end
if numel(iid)==1,    iid    = repmat(iid,1,n);    end
if numel(alpha)==1,  alpha  = repmat(alpha,1,n);  end

stat = zeros(1,n); pValue = zeros(1,n); cValue = zeros(1,n); ratio = zeros(1,n);
for i = 1:n
    q = period(i);

    % Lo-MacKinlay work with nq returns, so the sample is truncated to the
    % largest length whose return count is divisible by q. This matters only
    % when mod(T-1,q) ~= 0; for q dividing T-1 it is a no-op, which is why it
    % is easy to miss on power-of-two length data.
    nKeep = floor((T-1)/q)*q + 1;
    yq = y(1:nKeep);
    nq = nKeep - 1;
    muq = (yq(end)-yq(1))/nq;
    rq = diff(yq) - muq;
    saq = sum(rq.^2)/(nq-1);

    d = yq(q+1:end) - yq(1:end-q) - q*muq;
    m = q*(nq-q+1)*(1-q/nq);
    sc = sum(d.^2)/m;
    ratio(i) = sc/saq;

    if iid(i)
        theta = 2*(2*q-1)*(q-1)/(3*q*nq);
    else
        % Heteroskedasticity-consistent variance (Lo & MacKinlay 1988)
        den = sum(rq.^2)^2;
        theta = 0;
        for j = 1:q-1
            num = sum(rq(j+1:end).^2 .* rq(1:end-j).^2);
            theta = theta + (2*(q-j)/q)^2 * (num/den);
        end
    end
    stat(i) = (ratio(i)-1)/sqrt(theta);
    pValue(i) = 2*normcdf(-abs(stat(i)));   % avoids underflow in 1-normcdf
    cValue(i) = norminv(1-alpha(i)/2);
end
h = pValue < alpha;
end
