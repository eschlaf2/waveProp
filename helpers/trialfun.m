function [trl, event] = trialfun(cfg)

hdr = ft_read_header(cfg.dataset);
trl = round([cfg.starttime - cfg.padding(1), ...
	cfg.endtime + cfg.padding(2), ...
	cfg.padding(1)] * hdr.Fs);

event = [];
