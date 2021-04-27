function files = get_files(F)
N = height(F.Seizures);
fname =@(ii) sprintf('%s_Seizure%d_fits.mat', ...
    F.Seizures.patient{ii}, F.Seizures.seizure(ii));
files = arrayfun(@(ii) dir(fname(ii)), 1:N);
end
