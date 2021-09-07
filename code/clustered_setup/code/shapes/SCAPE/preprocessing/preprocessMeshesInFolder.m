files = dir('*.off');
for ii = 1:numel(files);
    preprocessMeshScape(files(ii).name, '', false, [])
end