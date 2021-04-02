function C = Func_map_to_bloc(dx, Nx, res, img)

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
for n = 1:N
    for m = 1:N
        mask = zeros(size(B));
        mask((n-1)*Nres+1:n*Nres, (m-1)*Nres+1:m*Nres) = 1;
        % imagesc(mask);
        avg = mean(nonzeros(B .* mask));
        B((n-1)*Nres+1:n*Nres, (m-1)*Nres+1:m*Nres) = avg;
        % pause(.2);
    end
end

C = B(Npad+1:N*Nres-Npad, Npad+1:N*Nres-Npad);

end