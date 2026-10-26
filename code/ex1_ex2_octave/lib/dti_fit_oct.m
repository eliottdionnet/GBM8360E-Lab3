function FitResults = dti_fit_oct(obj,data,nonlinear)
% Copie de dti.fit (qMRLab 3.0.0) ou l'etape non-lineaire (Rician) tourne aussi sous Octave.
if nargin<3, nonlinear = ~strcmp(obj.options.fittingtype,'linear'); end
Prot = ConvertSchemeUnits(obj.Prot.DiffusionData.Mat,0,1);
S0 = scd_preproc_getS0(data.DiffusionData,Prot);
if min(data.DiffusionData)<0, data.DiffusionData = max(0,data.DiffusionData); end
D = scd_model_dti(max(eps,data.DiffusionData)./max(eps,S0),Prot);
residue = 0;
if nonlinear
    if isfield(data,'SigmaNoise') && ~isempty(data.SigmaNoise)
        SigmaNoise = data.SigmaNoise(1);
    else
        SigmaNoise = computesigmanoise(obj.Prot.DiffusionData.Mat,data.DiffusionData);
    end
    if isempty(SigmaNoise), SigmaNoise = max(data.DiffusionData/100); end
    if SigmaNoise
        f = @(x) double(-2*sum(scd_model_likelihood_rician(data.DiffusionData,max(eps,S0.*equation(obj,x)),SigmaNoise)));
        [xopt,residue] = fminunc(f, D(:), optimset('MaxIter',20,'Display','off'));
        D(:) = xopt;
    end
end
[~,L] = eig(D); L = sort(diag(L),'descend'); L_mean = sum(L)/3;
FitResults.FA = sqrt(3/2)*sqrt(sum((L-L_mean).^2))/sqrt(sum(L.^2));
FitResults.L1 = L(1); FitResults.L2 = L(2); FitResults.L3 = L(3);
FitResults.D = D(:); FitResults.residue = residue; FitResults.S0 = S0(1);
end
