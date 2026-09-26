function y = normpdf(x, mu, sigma)
% Standard normal-density formula used by MATLAB's normpdf.
y = exp(-0.5 .* ((x-mu)./sigma).^2) ./ (sqrt(2*pi).*sigma);
end
