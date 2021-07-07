
function save_fr_dir(F, sub)
        
    assert(ismember(sub, 'BCDEH'));
    h = F.(sub);
    aa = h.Children;
    ylabel(aa, 'Direction');
    xlabel(aa, 'Time (s)');
    metric = h.Tag
    
    h_new = figure('position', [0 0 2.5 1.25]);
    aa_new = copyobj(aa, h_new);
    title(aa_new, F.MetricNames.(metric));
    aa_new.OuterPosition = [0 0 1 1];
    
    fname = [F.prefix 'bcde_' metric];
    F.print(h_new, fname);
    close(h_new);

end

