function C = Func_map_to_bloc_Exaustive(dx, Nx, res, img)

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

B = repmat(B, [1 1 Nres]);
for starting_pt = 1:Nres
    if starting_pt == 1
        for n = 1:N
            for m = 1:N
                mask = zeros(size(B));
                if m ~= N
                    mask((n-1)*Nres+1:n*Nres, ((m-1)*Nres+starting_pt):(m*Nres+starting_pt-1), starting_pt) = 1;
                else
                    mask((n-1)*Nres+1:n*Nres, ((m-1)*Nres+starting_pt):end, starting_pt) = 1;
                end
                avg = mean(nonzeros(B(:,:,starting_pt) .* mask(:,:,starting_pt)));

                if m ~= N
                    B((n-1)*Nres+1:n*Nres, ((m-1)*Nres+starting_pt):(m*Nres+starting_pt-1), starting_pt) = avg;
                else
                    B((n-1)*Nres+1:n*Nres, ((m-1)*Nres+starting_pt):end, starting_pt) = avg;
                end
                %if n ~= 1
                %    imagesc(mask(:,:,starting_pt));
                %    pause;
                %end
            end
        end
    else
        for n = 1:N
            for m = 1:N+1
                mask = zeros(size(B));
                if m ~= N+1 && m ~= 1
                    mask((n-1)*Nres+1:n*Nres, ((m-2)*Nres+starting_pt):((m-1)*Nres+starting_pt-1), starting_pt) = 1;
                elseif m == 1
                    mask((n-1)*Nres+1:n*Nres, ((m-1)*Nres+1):((m-1)*Nres+starting_pt-1), starting_pt) = 1;
                else
                    mask((n-1)*Nres+1:n*Nres, ((m-2)*Nres+starting_pt):end, starting_pt) = 1;
                end

                %size(B)
                %size(mask)
                %starting_pt
                %m
                avg = mean(nonzeros(B(:,:,starting_pt) .* mask(:,:,starting_pt)));
                if m ~= N+1 && m ~= 1
                    B((n-1)*Nres+1:n*Nres, ((m-2)*Nres+starting_pt):((m-1)*Nres+starting_pt-1), starting_pt) = avg;
                elseif m == 1
                    B((n-1)*Nres+1:n*Nres, ((m-1)*Nres+1):((m-1)*Nres+starting_pt-1), starting_pt) = avg;
                else
                    B((n-1)*Nres+1:n*Nres, ((m-2)*Nres+starting_pt):end, starting_pt) = avg;
                end

                % B(((n-1)*Nres+starting_pt):(n*Nres+starting_pt-1), ((m-1)*Nres+starting_pt):(m*Nres+starting_pt-1), starting_pt) = avg;
                % pause(.2);
                %if n ~= 1
                %    imagesc(mask(:,:,starting_pt));
                %    pause;
                %end
            end
        end
    end
end

C = B(Npad+1:N*Nres-Npad, Npad+1:N*Nres-Npad,:);

end