function Preview(sim)
    if isempty(sim.Qe_movie)
        load(sprintf('%s_%d_info.mat', sim.basename, sim.sim_num), 'Qe_movie');
        aspect_ratio = sim.grid_size(2) / sim.grid_size(1);
        for ii = 1:length(Qe_movie)
            Qe_movie(ii).cdata = ...
                imresize(Qe_movie(ii).cdata, [500 500 * aspect_ratio]);
        end
        sim.Qe_movie = Qe_movie;
    end
    h = figure(); fullwidth(true);
    T = tiledlayout(h, 10, 20);
    for ii = 1:min(numel(sim.Qe_movie), 200)
        nexttile(T, ii); 
        imagesc(sim.Qe_movie(ii).cdata); 
        colormap(sim.Qe_movie(ii).colormap);
        xticks([]); yticks([]);
        axis square;
    end
%     movie(h, sim.Qe_movie);
end