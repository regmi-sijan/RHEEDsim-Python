function close(varargin)
% Only accept the progress handle returned by this harness's waitbar.
if nargin ~= 1 || varargin{1} ~= -123456789
  error('Validation close shim received an unexpected handle');
end
end
