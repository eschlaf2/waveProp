function out = foldercontents2struct(path, filenames)
    % Load the contents of a directory into a structure where each
    % structure field is the name of a file
    % out = foldercontents2struct(path)
    if nargin < 2, filenames = ''; end
    
    files = {dir([path filesep filenames '*.mat']).name};
    files = cellfun(@(ss) strrep(ss, '.mat', ''), ...
        files, 'uni', 0);
    if isempty(files); out = []; end
    for ff = files
        field = load([path filesep ff{:}]);
        pn = fieldnames(field);
        out.(ff{:}) = field.(pn{:});
    end
    

end
