function [C, label_mat] = Func_map_to_bloc_3D(dx, dz, Nx, Nz, res, res_through, img)

Nres = round(res / dx); % res / dy
Nmod = mod(Nx,Nres);
if Nmod == 0
    Npad = 0;
    B = img;
    N = fix(Nx/Nres);
else
    Npad = (Nres - Nmod) / 2;
    disp(cat(2, 'Npad is: ', num2str(Npad)));
    N = fix(Nx/Nres) + 1;
    if fix(Npad) == Npad
        % temp = padarray(img, [Npad, Npad], 'symmetric', 'post');
        % B = padarray(temp, [Npad, Npad], 'symmetric', 'pre');
        B = padarray(img, [Npad, Npad], 'symmetric', 'both');
    else
        Npad1 = ceil(Npad);
        Npad2 = floor(Npad);
        temp = padarray(img, [Npad1, Npad1], 'symmetric', 'post');
        B = padarray(temp, [Npad2, Npad2], 'symmetric', 'pre');
    end
end

Nres_throughplane = round(res_through / dz); % res / dy
Nmod_z = mod(Nz, Nres_throughplane);
if Nmod_z == 0
    Npad_z = 0;
    B_z = B;
    N_3d = fix(Nz/Nres_throughplane);
else
    Npad_z = (Nres_throughplane - Nmod_z) / 2;
    %temp = padarray(B, [1, 1, Npad_z], 'symmetric', 'post');
    %B_z = padarray(temp, [1, 1, Npad_z], 'symmetric', 'pre');
    B_z = padarray(B, [1, 1, Npad_z], 'symmetric', 'both');
    N_3d = fix(Nz/Nres_throughplane) + 1;
end


% sliding windows
% figure();
sz = size(B_z);
label_mat = zeros(sz);
% B_z_temp = B_z;
count = 1;
for n = 1:N
    for m = 1:N
        %mask = zeros([sz(1), sz(2)]);
        %mask((n-1)*Nres+1:n*Nres, (m-1)*Nres+1:m*Nres) = 1;

        avg_xy = mean(B_z((n-1)*Nres+1:n*Nres, (m-1)*Nres+1:m*Nres,:), [1,2]); % 1x1xNz

        B_z_array = zeros(1, 1, sz(3));
        for l = 1:N_3d

            avg = mean(nonzeros(avg_xy(:, :, (l-1)*Nres_throughplane+1:l*Nres_throughplane)));
            % imagesc(mask);
            
            B_z_array(:,:,(l-1)*Nres_throughplane+1:l*Nres_throughplane) = avg;
            % pause(.2);

            % temp = ones(Nres, Nres, Nres_throughplane);
            label_mat((n-1)*Nres+1:n*Nres, (m-1)*Nres+1:m*Nres, (l-1)*Nres_throughplane+1:l*Nres_throughplane) = count;
            count = count + 1;
        end

        B_z((n-1)*Nres+1:n*Nres, (m-1)*Nres+1:m*Nres,:) = repmat(B_z_array, [Nres, Nres, 1]);
    end
end

if fix(Npad) == Npad
    C = B_z(Npad+1:N*Nres-Npad, Npad+1:N*Nres-Npad, Npad_z+1:N_3d*Nres_throughplane-Npad_z);
    label_mat = label_mat(Npad+1:N*Nres-Npad, Npad+1:N*Nres-Npad, Npad_z+1:N_3d*Nres_throughplane-Npad_z);
else
    C = B_z(Npad2+1:N*Nres-Npad1, Npad2+1:N*Nres-Npad1, Npad_z+1:N_3d*Nres_throughplane-Npad_z);
    label_mat = label_mat(Npad2+1:N*Nres-Npad1, Npad2+1:N*Nres-Npad1, Npad_z+1:N_3d*Nres_throughplane-Npad_z);
end

end