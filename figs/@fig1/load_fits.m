function load_fits(F)

metrics = struct2cell(F.Metrics)';
metrics(cellfun(@isempty, metrics)) = [];
fits = WaveProp.load(F.Files, metrics);

for ii = 1:length(F.Files)
    for mm = string(metrics)
        fits.(mm)(ii).MinFinite = 1;
        fits.(mm)(ii).Early = F.Early;
        fits.(mm)(ii).Late = F.Late;
        fits.(mm)(ii).RotateBy = F.true_angle(ii);
    end
end
F.Fits = fits;

end

