function C = Func_map_to_bloc_encoded2(dx, Nx, res, img)
% In order to find overlapping pixels with ground truth hemo map,
% The blocks are encoded
% So I only need encoded masks
% Extaustively slides the window
Nres = round(res / dx); % res / dy
Nmod = mod(Nx,Nres);
if Nmod == 0
    Npad = 0;
    B = img;
    N = fix(Nx/Nres);
else
    Npad = (Nres - mod(Nx,Nres)) / 2;
    temp = padarray(img, [Npad, Npad], 'symmetric', 'post');
    B = padarray(temp, [Npad, Npad], 'symmetric', 'pre');
    N = fix(Nx/Nres) + 1;
end

% sliding windows
% figure();
mask = zeros([size(B), Nres]);

for starting_pt = 1:Nres
    for n = 1:N
        for m = 1:N
            %mask = zeros(size(B));
            mask(((n-1)*Nres+starting_pt):(n*Nres+starting_pt-1), ((m-1)*Nres+starting_pt):(m*Nres+starting_pt-1), starting_pt) = N*(n-1)+m;
            % imagesc(mask);
            % avg = mean(nonzeros(B .* mask));
            % B((n-1)*Nres+1:n*Nres, (m-1)*Nres+1:m*Nres) = avg;
            % pause(.2);
        end
    end
end

C = mask(Npad+1:N*Nres-Npad, Npad+1:N*Nres-Npad, :);

end