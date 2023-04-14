function C = Func_map_to_bloc_encoded_exaustive(dx, Nx, res, img)
% In order to find overlapping pixels with ground truth hemo map,
% The blocks are encoded
% So I only need encoded masks
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
mask = zeros(size(B));

for starting_pt = 1:Nres
    if starting_pt == 1
        for n = 1:N
            for m = 1:N
                if m ~= N
                    mask((n-1)*Nres+1:n*Nres, ((m-1)*Nres+starting_pt):(m*Nres+starting_pt-1), starting_pt) = N*(n-1)+m;
                else
                    mask((n-1)*Nres+1:n*Nres, ((m-1)*Nres+starting_pt):end, starting_pt) = N*(n-1)+m;
                end
            end
        end
    else
        for n = 1:N
            for m = 1:N+1
                if m ~= N+1 && m ~= 1
                    mask((n-1)*Nres+1:n*Nres, ((m-2)*Nres+starting_pt):((m-1)*Nres+starting_pt-1), starting_pt) = N*(n-1)+m;
                elseif m == 1
                    mask((n-1)*Nres+1:n*Nres, ((m-1)*Nres+1):((m-1)*Nres+starting_pt-1), starting_pt) = N*(n-1)+m;
                else
                    mask((n-1)*Nres+1:n*Nres, ((m-2)*Nres+starting_pt):end, starting_pt) = N*(n-1)+m;
                end
            end
        end
    end
end

C = mask(Npad+1:N*Nres-Npad, Npad+1:N*Nres-Npad,:);


end