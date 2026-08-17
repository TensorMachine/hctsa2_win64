function d = BF_iddata(y,u,Ts)
% BF_iddata   Toolbox-free stand-in for iddata().
%
% hctsa only ever calls iddata(y,[],1): a single-output time series with no
% input and unit sample time. That carries no information beyond the vector
% itself, so this returns a plain column vector.
%
% Returning a vector rather than a struct is deliberate. The MF_* operations
% slice these objects (`y(r(i,1):r(i,2))` in MF_FitSubsegments), which a scalar
% struct cannot support, and they read `.y` off them, which a vector cannot.
% A vector keeps slicing, length(), size() and arithmetic all working; the `.y`
% reads were changed to plain identifiers instead. The alternative -- a classdef
% with subsref overloading -- would avoid the caller edits but adds a class to a
% codebase that is otherwise entirely functions and structs.

if nargin >= 2 && ~isempty(u)
    error('BF_iddata:inputNotSupported', ...
          'Input signals are not supported: hctsa only uses output-only data.');
end
if nargin >= 3 && ~isempty(Ts) && Ts ~= 1
    error('BF_iddata:sampleTime','Only Ts = 1 is supported.');
end
d = y(:);
end
