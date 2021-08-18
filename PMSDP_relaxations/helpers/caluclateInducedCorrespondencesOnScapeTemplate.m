function corrToSave = caluclateInducedCorrespondencesOnScapeTemplate(first,second,m1idxfull,m2idxfull)

% load template mesh
[V1 F1] = read_off(sprintf('mesh%03d.off',first));
[V2 F2] = read_off(sprintf('mesh%03d.off',second));
le=reshape(read(sprintf('test_%03d.txt',first)),3,[]);
ri=reshape(read(sprintf('test_%03d.txt',second)),3,[]);
% find correspondences
fullScapeM1Idx = knnsearch(V1',le(:,m1idxfull)');
fullScapeM2Idx = knnsearch(V2',ri(:,m2idxfull)');
corrToSave = [fullScapeM1Idx-1 fullScapeM2Idx-1]';
corrToSave = corrToSave(:);

end