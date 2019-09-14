function name = checkname(name)

strmap = {...
    '.', 'p'; ...
    '-', 'M'; ...
    };

for ii = 1:size(strmap, 1)
    name = strrep(name, strmap{ii, 1}, strmap{ii, 2}); 
end
