function p = pdf(name, x, s, sigma)
% Remplacement minimal de pdf('rician',...) (Statistics Toolbox) pour Octave.
if ~strcmpi(name,'rician'), error('pdf: seul rician est implemente'); end
sigsq = sigma.^2;
p = (x./sigsq) .* exp(-0.5*(x-s).^2./sigsq) .* besseli(0, x.*s./sigsq, 1);
p(x<0) = 0;
end
