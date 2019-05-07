function [] = mea_time_to_function(mea)

Time = mea.Time;
temp = @() linspace(Time(1), Time(end), numel(Time));
mea.Time = temp;