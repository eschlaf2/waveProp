function Preview(sim)
    if isempty(sim.Qe_movie)
        load(sprintf('%s_%d_info.mat', sim.basename, sim.sim_num), 'Qe_movie', 'K_movie');
        aspect_ratio = sim.grid_size(2) / sim.grid_size(1);
        for ii = 1:length(Qe_movie)
            Qe_movie(ii).cdata = ...
                imresize(Qe_movie(ii).cdata, [500 500 * aspect_ratio]);
            K_movie(ii).cdata = ...
                imresize(K_movie(ii).cdata, [500 500 * aspect_ratio]);
        end
        sim.Qe_movie = Qe_movie;
        sim.K_movie = K_movie;
    end

     h = figure('units', 'inches', ...
        'position', [0    0.8472   10.6250   10.2083], ...
        'name', 'K_movie'); 
    T = tiledlayout(h, 10, 10, 'tilespacing', 'compact');
    for ii = 1:min(numel(sim.Qe_movie), 100)
        nexttile(T, ii); 
        imagesc(sim.Qe_movie(ii).cdata); 
        colormap(sim.Qe_movie(ii).colormap);
        xticks([]); yticks([]);
        axis square;
    end
    
    h = figure('units', 'inches', ...
        'position', [0    0.8472   10.6250   10.2083], ...
        'name', 'K_movie'); 
    T = tiledlayout(h, 10, 10, 'tilespacing', 'compact');
    for ii = 1:min(numel(sim.K_movie), 100)
        nexttile(T, ii); 
        image(sim.K_movie(ii).cdata); 
        colormap(sim.K_movie(ii).colormap);
        xticks([]); yticks([]);
        axis square;
    end
    
%     movie(h, sim.Qe_movie);
end