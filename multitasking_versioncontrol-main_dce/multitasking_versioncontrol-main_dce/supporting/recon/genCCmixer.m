function [mixer, eigValues] = genCCmixer(roi,roni)

[eigVector, eigValues] = eig(roi,roni);
[~, idx]  = sort(sqrt(abs(diag(eigValues))),'descend');
mixer = eigVector(:,idx);

% [eigVector, eigValues] = eig(roni\roi);
% [~, idx] = sort(sqrt(abs(diag(eigValues))),'descend');
% mixer = eigVector(:,idx);

end