function [water, fat1, fat2, freq, R2s, iter, model, fitting_error_matrix] = fit_MP_IDEAL_R2_Chris(s0, t, f_fat, f0, R2s, max_iter, Mask)
matrix_size = size(s0);
numvox = prod(matrix_size(1:end-1));
numte = matrix_size(end);

if nargin<6
    max_iter = 30;
end
if nargin <5
    R2s = zeros([1 numvox]);
end
if nargin<4
    f0 = zeros([1 numvox]);
end
if nargin<3
    f_fat = -420;
end

if numel(f0) == 1
    f0 = f0*ones([1 numvox]);
end
if numel(R2s) == 1
    R2s = R2s*ones([1 numvox]);
end

s0 = permute(s0, [length(matrix_size) 1:length(matrix_size)-1]);
s0 = reshape(s0,[numte numvox]);
R2s = reshape(R2s,[1 numvox]);
f0 = reshape(f0,[1 numvox]);
t = reshape(t,[numte 1]);
t = repmat(t,[1 numvox]);

% O = ones([numte numvox]).*exp(repmat(-R2s,[numte 1]).*t);
% C = real(exp(-1i*2*pi*f_fat*t)).*exp(repmat(-R2s, [numte 1]).*t);

O = ones([numte numvox]);
C1 = exp(-1i*2*pi*f_fat(1)*t);
C2 = exp(-1i*2*pi*f_fat(2)*t);
% C3 = exp(-1i*2*pi*f_fat(3)*t);

% y = zeros([5 numvox]);
y = zeros([4 numvox]);
y(1,:) = 1;
% dy = zeros([5 numvox]);
dy = zeros([4 numvox]);
dy(1,:) = 1e4;
iter = 0;

y(1,:) = f0+1i*R2s/(2*pi);

P = exp(-1i*2*pi*repmat(y(1,:),[numte 1]).*t);  %complex phasor
% y(2:5,:) = invA4(P.*O, P.*C1, P.*C2, P.*C3, s0);
y(2:4,:) = invA3(P.*O, P.*C1, P.*C2, s0);
% dy_momentum = zeros(size(y));
% beta = 0.9;
update = dy(1,:);
while (iter<max_iter)&&(sqrt(sum(real(update).^2,2)/numvox)>0.1)
    % sn = P.*O.*repmat(y(2,:),[numte 1]) + P.*C1.*repmat(y(3,:),[numte 1]) + P.*C2.*repmat(y(4,:),[numte 1]) + P.*C3.*repmat(y(5,:),[numte 1]);
    sn = P.*O.*repmat(y(2,:),[numte 1]) + P.*C1.*repmat(y(3,:),[numte 1]) + P.*C2.*repmat(y(4,:),[numte 1]);
    sr = s0 - sn;
%     gr = -2*pi*t.*(-repmat(y(2,:),[numte 1]).*z - repmat(y(3,:),[numte 1]).*O - repmat(y(4,:),[numte 1]).*d - repmat(y(5,:),[numte 1]).*C);
%     gi = -2*pi*t.*(repmat(y(2,:),[numte 1]).*O - repmat(y(3,:),[numte 1]).*z + repmat(y(4,:),[numte 1]).*C - repmat(y(5,:),[numte 1]).*d);

    Bcol01 = -1i*2*pi*t.*sn;
    dy = invA4(Bcol01, P.*O, P.*C1, P.*C2, sr); 
    % dy = invA5(Bcol01, P.*O, P.*C1, P.*C2, P.*C3, sr);
    y = y+dy;
%     dy = invA4(Bcol01, P.*O, P.*C1, P.*C2, sr);    
%     dy_momentum = (1-beta).*dy+0.0*dy_momentum;
%     y = y+dy_momentum;
    iter = iter+1;

%     temp = y(1,:);
%     temp(abs(imag(temp))*2*pi>50)=real(temp(abs(imag(temp))*2*pi>50));
%     y(1,:) = temp;  %hann_low
    
    P = exp(-1i*2*pi*repmat(y(1,:),[numte 1]).*t);  %complex phasor
    % y(2:5,:) = invA4(P.*O, P.*C1, P.*C2, P.*C3, s0);
    y(2:4,:) = invA3(P.*O, P.*C1, P.*C2, s0);

    update = dy(1,:);
    update(isnan(update)) = 0;
    update(isinf(update)) = 0;
%     error_iter(iter) = sqrt(sum(abs(update(:)).^2)/numvox);

    model_temp = permute(sn, [2 1]);
    numvoxnnz = nnz(Mask);
    fitting_error_vector= sqrt(sum((abs(model_temp) - abs(s0')).^2,2))./sum(abs(s0'),2);
    fitting_error_vector(isnan(fitting_error_vector)) = 0;
%     numvoxnnz = 1;%nnz(Mask);
    error_iter(iter) = sqrt(sum(fitting_error_vector(Mask~=0).^2,1))/numvoxnnz;
%     error_iter(iter) = sum(abs(sr(Mask~=0)));
end

figure;plot(error_iter);xlim([0,iter]);%ylim([0,100]);

freq = reshape(real(y(1,:)),matrix_size(1:end-1));
R2s = -reshape(imag(y(1,:))*2*pi,matrix_size(1:end-1));
water = reshape(y(2,:),matrix_size(1:end-1));
fat1 = reshape(y(3,:),matrix_size(1:end-1));
fat2 = reshape(y(4,:),matrix_size(1:end-1));
% fat3 = reshape(y(5,:),matrix_size(1:end-1));

% model = P.*O.*repmat(y(2,:),[numte 1]) + P.*C1.*repmat(y(3,:),[numte 1])+ P.*C2.*repmat(y(4,:),[numte 1])+ P.*C3.*repmat(y(5,:),[numte 1]);
model = P.*O.*repmat(y(2,:),[numte 1]) + P.*C1.*repmat(y(3,:),[numte 1])+ P.*C2.*repmat(y(4,:),[numte 1]);
% model = reshape(model,matrix_size);
model = permute(model, [2 1]);
% numvoxnnz = nnz(Mask);

fitting_error_vector= sqrt(sum((abs(model) - abs(s0')).^2,2))./sum(abs(s0'),2);

numvoxnnz = 1;%nnz(Mask);
fitting_error = sqrt(sum(fitting_error_vector(Mask~=0).^2,1))/numvoxnnz;
fitting_error_matrix = reshape(fitting_error_vector,matrix_size(1:end-1));

model = reshape(model,matrix_size);
% s_model = [s0;model];
% figure; plot(s0,'ro'); hold on; plot(model,'bx');hold off;
% axis([-max(abs(real(s_model))) max(abs(real(s_model))) -max(abs(imag(s_model))) max(abs(imag(s_model)))]*1.2)

freq(isinf(freq)) = 0;
freq(isnan(freq)) = 0;
freq(abs(freq)>10e4) = 0;
R2s(isinf(R2s)) = 0;
R2s(isnan(R2s)) = 0;
water(isinf(water)) = 0;
water(isnan(water)) = 0;
water(abs(water)>10e5) = 0;
fat1(isinf(fat1)) = 0;
fat1(isnan(fat1)) = 0;
fat1(abs(fat1)>10e5) = 0;
fat2(isinf(fat2)) = 0;
fat2(isnan(fat2)) = 0;
fat2(abs(fat2)>10e5) = 0;
%fat3(isinf(fat3)) = 0;
%fat3(isnan(fat3)) = 0;
%fat3(abs(fat3)>10e5) = 0;
model(isinf(model)) = 0;
model(isnan(model)) = 0;

% freq(isinf(freq)) = 0;
% freq(isnan(freq)) = 0;
% 
% R2s(isinf(R2s)) = 0;
% R2s(isnan(R2s)) = 0;
% water(isinf(water)) = 0;
% water(isnan(water)) = 0;
% 
% fat(isinf(fat)) = 0;
% fat(isnan(fat)) = 0;
% 
% model(isinf(model)) = 0;
% model(isnan(model)) = 0;
end

function x=invA2(col1, col2, y)
% assemble A^H*A
a11 = sum(conj(col1).*col1, 1);
a12 = sum(conj(col1).*col2, 1);
a22 = sum(conj(col2).*col2, 1);

% inversion of A^H*A
d = (a11.*a22 - a12.*conj(a12));
ia11 = a22./d;
ia12 = -a12./d;
ia22 = a11./d;

% y project onto A^H
py1 = sum(conj(col1).*y,1);
py2 = sum(conj(col2).*y,1);
% calculate x
x(1,:) = sum(ia11.*py1 + ia12.*py2, 1);
x(2,:) = sum(conj(ia12).*py1 + ia22.*py2, 1);
end

function x=invA3(col1, col2, col3, y)
% assemble B^H*B
b11 = sum(conj(col1).*col1,1);
b12 = sum(conj(col1).*col2,1);
b13 = sum(conj(col1).*col3,1);
b22 = sum(conj(col2).*col2,1);
b23 = sum(conj(col2).*col3,1);
b33 = sum(conj(col3).*col3,1);

% inversion of B'*B
d = (b13.*conj(b12).*conj(b23) + b11.*b22.*b33 + b12.*b23.*conj(b13) - b13.*b22.*conj(b13) - b11.*b23.*conj(b23) - b12.*b33.*conj(b12));
ib11 = (b22.*b33 - b23.*conj(b23))./d;
ib12 = -(b12.*b33 - b13.*conj(b23))./d;
ib13 = (b12.*b23 - b13.*b22)./d;
ib22 = (b11.*b33 - b13.*conj(b13))./d;
ib23 = -(b11.*b23 - b13.*conj(b12))./d;
ib33 =  (b11.*b22 - b12.*conj(b12))./d;

% y project onto B'
py1 = sum(conj(col1).*y,1);
py2 = sum(conj(col2).*y,1);
py3 = sum(conj(col3).*y,1);
% calculate x
x(1,:) = sum(ib11.*py1 + ib12.*py2 + ib13.*py3, 1);
x(2,:) = sum(conj(ib12).*py1 + ib22.*py2 + ib23.*py3, 1);
x(3,:) = sum(conj(ib13).*py1 + conj(ib23).*py2 + ib33.*py3, 1);
end

function x=invA4(col1, col2, col3, col4, y)
% assemble A'*A
a11 = sum(conj(col1).*col1,1);
a12 = sum(conj(col1).*col2,1);
a13 = sum(conj(col1).*col3,1);
a14 = sum(conj(col1).*col4,1);
a22 = sum(conj(col2).*col2,1);
a23 = sum(conj(col2).*col3,1);
a24 = sum(conj(col2).*col4,1);
a33 = sum(conj(col3).*col3,1);
a34 = sum(conj(col3).*col4,1);
a44 = sum(conj(col4).*col4,1);

% inversion of A'*A
%{
% d = (a33.*a44.*a12.^2 - a12.^2.*a34.^2 - 2.*a44.*a12.*a13.*a23 + 2.*a12.*a13.*a24.*a34 + 2.*a12.*a14.*a23.*a34 - 2.*a33.*a12.*a14.*a24 - a13.^2.*a24.^2 + a22.*a44.*a13.^2 + 2.*a13.*a14.*a23.*a24 - 2.*a22.*a13.*a14.*a34 - a14.^2.*a23.^2 + a22.*a33.*a14.^2 + a11.*a44.*a23.^2 - 2.*a11.*a23.*a24.*a34 + a11.*a33.*a24.^2 + a11.*a22.*a34.^2 - a11.*a22.*a33.*a44);
% ia11 = (a44.*a23.^2 - 2.*a23.*a24.*a34 + a33.*a24.^2 + a22.*a34.^2 - a22.*a33.*a44)./d;
% ia12 = -(a12.*a34.^2 - a13.*a24.*a34 - a14.*a23.*a34 + a14.*a24.*a33 + a13.*a23.*a44 - a12.*a33.*a44)./d;
% ia13 = -(a13.*a24.^2 - a14.*a23.*a24 - a12.*a24.*a34 + a14.*a22.*a34 + a12.*a23.*a44 - a13.*a22.*a44)./d;
% ia14 = -(a14.*a23.^2 - a13.*a23.*a24 - a12.*a23.*a34 + a12.*a24.*a33 + a13.*a22.*a34 - a14.*a22.*a33)./d;
% ia22 = (a44.*a13.^2 - 2.*a13.*a14.*a34 + a33.*a14.^2 + a11.*a34.^2 - a11.*a33.*a44)./d;
% ia23 = -(a14.^2.*a23 - a13.*a14.*a24 - a12.*a14.*a34 + a11.*a24.*a34 + a12.*a13.*a44 - a11.*a23.*a44)./d;
% ia24 = -(a13.^2.*a24 - a13.*a14.*a23 - a12.*a13.*a34 + a12.*a14.*a33 + a11.*a23.*a34 - a11.*a24.*a33)./d;
% ia33 = (a44.*a12.^2 - 2.*a12.*a14.*a24 + a22.*a14.^2 + a11.*a24.^2 - a11.*a22.*a44)./d;
% ia34 = -(a12.^2.*a34 - a12.*a13.*a24 - a12.*a14.*a23 + a13.*a14.*a22 + a11.*a23.*a24 - a11.*a22.*a34)./d;
% ia44 = (a33.*a12.^2 - 2.*a12.*a13.*a23 + a22.*a13.^2 + a11.*a23.^2 - a11.*a22.*a33)./d;
%}
d = (a11.*a22.*a33.*a44 + a11.*a23.*a34.*conj(a24) + a11.*a24.*conj(a23).*conj(a34)...
     -a11.*a24.*a33.*conj(a24) - a11.*a23.*conj(a23).*a44 - a11.*a22.*a34.*conj(a34)...
     -a12.*conj(a12).*a33.*a44 - a13.*conj(a12).*a34.*conj(a24) - a14.*conj(a12).*conj(a23).*conj(a34)...
     +a14.*conj(a12).*a33.*conj(a24) + a13.*conj(a12).*conj(a23).*a44 + a12.*conj(a12).*a34.*conj(a34)...
     +a12.*a23.*conj(a13).*a44 + a13.*a24.*conj(a13).*conj(a24) + a14.*a22.*conj(a13).*conj(a34)...
     -a14.*a23.*conj(a13).*conj(a24) - a13.*a22.*conj(a13).*a44 - a12.*a24.*conj(a13).*conj(a34)...
     -a12.*a23.*a34.*conj(a14) - a13.*a24.*conj(a23).*conj(a14) - a14.*a22.*a33.*conj(a14)...
     +a14.*a23.*conj(a23).*conj(a14) + a13.*a22.*a34.*conj(a14) + a12.*a24.*a33.*conj(a14));

ia11 = (a22.*a33.*a44 +  a23.*a34.*conj(a24) + a24.*conj(a23).*conj(a34)...
        - a24.*a33.*conj(a24) - a23.*conj(a23).*a44 - a22.*a34.*conj(a34))./d;
ia12 = -(a12.*a33.*a44 +  a13.*a34.*conj(a24) + a14.*conj(a23).*conj(a34)...
        - a14.*a33.*conj(a24) - a13.*conj(a23).*a44 - a12.*a34.*conj(a34))./d;
ia13 = (a12.*a23.*a44 +  a13.*a24.*conj(a24) + a14.*a22.*conj(a34)...
        - a14.*a23.*conj(a24) - a13.*a22.*a44 - a12.*a24.*conj(a34))./d;
ia14 = -(a12.*a23.*a34 +  a13.*a24.*conj(a23) + a14.*a22.*a33...
        - a14.*a23.*conj(a23) - a13.*a22.*a34 - a12.*a24.*a33)./d;
ia22 = (a11.*a33.*a44 +  a13.*a34.*conj(a14) + a14.*conj(a13).*conj(a34)...
        - a14.*a33.*conj(a14) - a13.*conj(a13).*a44 - a11.*a34.*conj(a34))./d;
ia23 = -(a11.*a23.*a44 +  a13.*a24.*conj(a14) + a14.*conj(a12).*conj(a34)...
        - a14.*a23.*conj(a14) - a13.*conj(a12).*a44 - a11.*a24.*conj(a34))./d;
ia24 = (a11.*a23.*a34 +  a13.*a24.*conj(a13) + a14.*conj(a12).*a33...
        - a14.*a23.*conj(a13) - a13.*conj(a12).*a34 - a11.*a24.*a33)./d;
ia33 = (a11.*a22.*a44 +  a12.*a24.*conj(a14) + a14.*conj(a12).*conj(a24)...
        - a14.*a22.*conj(a14) - a12.*conj(a12).*a44 - a11.*a24.*conj(a24))./d;
ia34 = -(a11.*a22.*a34 +  a12.*a24.*conj(a13) + a14.*conj(a12).*conj(a23)...
        - a14.*a22.*conj(a13) - a12.*conj(a12).*a34 - a11.*a24.*conj(a23))./d;
ia44 = (a11.*a22.*a33 +  a12.*a23.*conj(a13) + a13.*conj(a12).*conj(a23)...
        - a13.*a22.*conj(a13) - a12.*conj(a12).*a33 - a11.*a23.*conj(a23))./d;
% y project onto A'
py1 = sum(conj(col1).*y,1);
py2 = sum(conj(col2).*y,1);
py3 = sum(conj(col3).*y,1);
py4 = sum(conj(col4).*y,1);
% calculate x
x(1,:) = sum(ia11.*py1 + ia12.*py2 + ia13.*py3 + ia14.*py4,1);
x(2,:) = sum(conj(ia12).*py1 + ia22.*py2 + ia23.*py3 + ia24.*py4,1);
x(3,:) = sum(conj(ia13).*py1 + conj(ia23).*py2 + ia33.*py3 + ia34.*py4,1);
x(4,:) = sum(conj(ia14).*py1 + conj(ia24).*py2 + conj(ia34).*py3 + ia44.*py4,1);
end

function x=invA5(col1, col2, col3, col4, col5,y)
% assemble B'*B
b11 = sum(conj(col1).*col1,1);
b12 = sum(conj(col1).*col2,1);b21 = conj(b12);
b13 = sum(conj(col1).*col3,1);b31 = conj(b13);
b14 = sum(conj(col1).*col4,1);b41 = conj(b14);
b15 = sum(conj(col1).*col5,1);b51 = conj(b15);
b22 = sum(conj(col2).*col2,1);
b23 = sum(conj(col2).*col3,1);b32 = conj(b23);
b24 = sum(conj(col2).*col4,1);b42 = conj(b24);
b25 = sum(conj(col2).*col5,1);b52 = conj(b25);
b33 = sum(conj(col3).*col3,1);
b34 = sum(conj(col3).*col4,1);b43 = conj(b34);
b35 = sum(conj(col3).*col5,1);b53 = conj(b35);
b44 = sum(conj(col4).*col4,1);
b45 = sum(conj(col4).*col5,1);b54 = conj(b45);
b55 = sum(conj(col5).*col5,1);
%{
% inversion of B'*B
d = (-b23.*b45.*b12.*b13.*b45-b11.*b23.*b45.*b34.*b25+b23.*b35.*b12.*b14.*b45+b23.*b45.*b14.*b13.*b25-b23.*b11.*b45.*b24.*b35+b11.*b23.*b44.*b35.*b25-b23.*b34.*b12.*b14.*b55+b23.*b44.*b12.*b13.*b55-b23.*b14.*b13.*b24.*b55+b23.*b15.*b14.*b34.*b25+b23.*b15.*b13.*b24.*b45-b23.*b44.*b35.*b12.*b15-b23.*b34.*b24.*b15.*b15+b23.*b45.*b34.*b12.*b15+b23.*b35.*b24.*b14.*b15-b23.*b44.*b15.*b13.*b25-b23.*b14.*b14.*b35.*b25-b11.*b23.*b35.*b24.*b45-b45.*b22.*b34.*b13.*b15-b45.*b33.*b14.*b12.*b25+b35.*b24.*b13.*b12.*b45-b35.*b24.*b14.*b12.*b35-b22.*b35.*b34.*b14.*b15+b22.*b33.*b45.*b14.*b15+b45.*b34.*b13.*b12.*b25-b45.*b33.*b24.*b12.*b15-b45.*b22.*b14.*b13.*b35+b45.*b22.*b13.*b13.*b45-b33.*b25.*b24.*b14.*b15-b22.*b33.*b44.*b15.*b15+b33.*b24.*b24.*b15.*b15+b45.*b33.*b12.*b12.*b45-b35.*b24.*b24.*b13.*b15+b45.*b24.*b12.*b13.*b35-b45.*b34.*b12.*b12.*b35+b34.*b24.*b15.*b12.*b35-b44.*b35.*b13.*b12.*b25-b44.*b25.*b12.*b13.*b35+b44.*b33.*b25.*b12.*b15-b45.*b13.*b13.*b24.*b25+b22.*b34.*b34.*b15.*b15+b34.*b24.*b25.*b13.*b15-b34.*b24.*b13.*b12.*b55+b44.*b35.*b12.*b12.*b35+b44.*b22.*b15.*b13.*b35+b34.*b12.*b12.*b34.*b55+b22.*b34.*b13.*b14.*b55+b22.*b14.*b13.*b34.*b55+b33.*b14.*b12.*b24.*b55+b44.*b13.*b13.*b25.*b25-b44.*b33.*b12.*b12.*b55+b24.*b34.*b15.*b13.*b25-b24.*b34.*b12.*b13.*b55+b24.*b34.*b35.*b12.*b15+b44.*b22.*b35.*b13.*b15+b14.*b12.*b34.*b35.*b25+b33.*b24.*b12.*b14.*b55-b25.*b34.*b14.*b13.*b25+b25.*b34.*b12.*b13.*b45-b15.*b14.*b22.*b34.*b35-b15.*b14.*b33.*b24.*b25-b33.*b25.*b12.*b14.*b45+b14.*b14.*b22.*b35.*b35+b14.*b14.*b33.*b25.*b25+b44.*b33.*b15.*b12.*b25-b44.*b22.*b13.*b13.*b55-b22.*b33.*b14.*b14.*b55+b14.*b13.*b24.*b25.*b35-b13.*b13.*b25.*b24.*b45-b25.*b13.*b14.*b34.*b25-b22.*b35.*b13.*b14.*b45-b15.*b12.*b34.*b34.*b25-b35.*b12.*b12.*b34.*b45-b35.*b12.*b14.*b24.*b35-b15.*b13.*b24.*b24.*b35-b33.*b15.*b12.*b24.*b45-b25.*b34.*b34.*b12.*b15-b22.*b15.*b13.*b34.*b45+b22.*b33.*b15.*b14.*b45+b13.*b13.*b24.*b24.*b55+b24.*b13.*b14.*b35.*b25+b34.*b12.*b14.*b25.*b35+b11.*b23.*b34.*b24.*b55+b11.*b22.*b33.*b44.*b55-b11.*b34.*b24.*b25.*b35+b11.*b45.*b33.*b24.*b25+b11.*b25.*b34.*b34.*b25+b11.*b45.*b22.*b34.*b35+b11.*b33.*b25.*b24.*b45-b11.*b44.*b22.*b35.*b35-b11.*b33.*b24.*b24.*b55+b11.*b22.*b35.*b34.*b45+b11.*b35.*b24.*b24.*b35-b11.*b22.*b34.*b34.*b55-b11.*b44.*b33.*b25.*b25-b11.*b22.*b33.*b45.*b45-b11.*b24.*b34.*b35.*b25-b23.*b11.*b23.*b44.*b55+b23.*b11.*b23.*b45.*b45+b23.*b11.*b44.*b25.*b35-b23.*b11.*b25.*b34.*b45+b23.*b11.*b24.*b34.*b55+b23.*b23.*b14.*b14.*b55-b23.*b23.*b15.*b14.*b45-b23.*b23.*b45.*b14.*b15+b23.*b23.*b44.*b15.*b15-b23.*b44.*b15.*b12.*b35+b23.*b15.*b12.*b34.*b45-b23.*b14.*b12.*b34.*b55+b23.*b25.*b34.*b14.*b15+b23.*b25.*b13.*b14.*b45+b23.*b15.*b14.*b24.*b35-b23.*b44.*b25.*b13.*b15-b23.*b24.*b13.*b14.*b55-b23.*b14.*b14.*b25.*b35-b23.*b45.*b13.*b12.*b45+b23.*b45.*b14.*b12.*b35+b23.*b44.*b13.*b12.*b55+b23.*b45.*b24.*b13.*b15-b23.*b24.*b34.*b15.*b15);
ib11=(-b23.*b23.*b44.*b55+b23.*b23.*b45.*b45+b23.*b44.*b25.*b35-b23.*b25.*b34.*b45-b23.*b45.*b24.*b35+b23.*b24.*b34.*b55-b23.*b45.*b34.*b25-b23.*b35.*b24.*b45+b23.*b44.*b35.*b25+b23.*b34.*b24.*b55+b22.*b33.*b44.*b55-b34.*b24.*b25.*b35+b45.*b33.*b24.*b25+b25.*b34.*b34.*b25+b45.*b22.*b34.*b35+b33.*b25.*b24.*b45-b44.*b22.*b35.*b35-b33.*b24.*b24.*b55+b22.*b35.*b34.*b45+b35.*b24.*b24.*b35-b22.*b34.*b34.*b55-b44.*b33.*b25.*b25-b22.*b33.*b45.*b45-b24.*b34.*b35.*b25)./d;
ib12=(-b14.*b24.*b35.*b35+b14.*b24.*b55.*b33-b14.*b55.*b23.*b34-b14.*b25.*b45.*b33+b14.*b35.*b25.*b34+b14.*b23.*b45.*b35+b24.*b35.*b45.*b13+b24.*b34.*b15.*b35-b24.*b45.*b15.*b33-b24.*b55.*b34.*b13-b12.*b34.*b45.*b35-b35.*b44.*b13.*b25+b45.*b12.*b45.*b33+b44.*b12.*b35.*b35+b55.*b12.*b34.*b34-b34.*b15.*b25.*b34+b55.*b23.*b44.*b13+b25.*b45.*b34.*b13-b23.*b44.*b15.*b35-b23.*b45.*b45.*b13-b35.*b45.*b34.*b12+b45.*b15.*b23.*b34-b55.*b44.*b12.*b33+b44.*b25.*b15.*b33)./d;
ib13=(-b14.*b23.*b55.*b24+b14.*b23.*b25.*b45+b14.*b35.*b25.*b24-b14.*b25.*b34.*b25-b14.*b22.*b35.*b45+b14.*b22.*b34.*b55+b23.*b45.*b15.*b24+b23.*b55.*b44.*b12-b23.*b44.*b25.*b15-b23.*b45.*b12.*b45+b24.*b13.*b24.*b55-b24.*b45.*b13.*b25-b24.*b35.*b24.*b15-b13.*b25.*b24.*b45+b35.*b24.*b12.*b45-b44.*b22.*b13.*b55+b44.*b22.*b35.*b15+b45.*b22.*b13.*b45-b44.*b35.*b12.*b25+b34.*b24.*b25.*b15+b45.*b34.*b12.*b25+b44.*b13.*b25.*b25-b45.*b22.*b34.*b15-b34.*b24.*b12.*b55)./d;
ib14=-(b14.*b23.*b35.*b25-b14.*b23.*b55.*b23+b14.*b55.*b22.*b33-b14.*b33.*b25.*b25+b14.*b35.*b25.*b23-b14.*b35.*b22.*b35-b23.*b34.*b15.*b25-b23.*b35.*b12.*b45+b23.*b55.*b12.*b34+b23.*b45.*b15.*b23-b45.*b15.*b22.*b33-b55.*b33.*b12.*b24-b55.*b22.*b13.*b34+b34.*b15.*b22.*b35-b35.*b24.*b15.*b23-b45.*b13.*b25.*b23+b55.*b13.*b23.*b24-b35.*b25.*b12.*b34-b35.*b13.*b25.*b24+b34.*b13.*b25.*b25+b35.*b35.*b12.*b24+b33.*b25.*b12.*b45+b35.*b22.*b13.*b45+b33.*b24.*b15.*b25)./d;
ib15=-(b23.*b23.*b14.*b45-b23.*b23.*b44.*b15+b23.*b44.*b12.*b35+b23.*b24.*b34.*b15-b23.*b12.*b34.*b45-b23.*b14.*b24.*b35-b23.*b13.*b24.*b45+b23.*b34.*b24.*b15-b23.*b14.*b34.*b25+b23.*b44.*b13.*b25-b22.*b34.*b34.*b15-b33.*b24.*b24.*b15-b22.*b33.*b14.*b45-b44.*b33.*b12.*b25+b22.*b13.*b34.*b45+b22.*b33.*b44.*b15+b14.*b22.*b34.*b35+b33.*b12.*b24.*b45-b24.*b34.*b13.*b25-b44.*b22.*b13.*b35+b14.*b33.*b24.*b25+b13.*b24.*b24.*b35-b34.*b24.*b12.*b35+b12.*b34.*b34.*b25)./d;
ib22=(-b14.*b13.*b45.*b35+b14.*b13.*b55.*b34+b14.*b45.*b15.*b33+b14.*b14.*b35.*b35-b14.*b55.*b33.*b14-b14.*b15.*b35.*b34-b13.*b45.*b15.*b34+b13.*b44.*b15.*b35+b13.*b45.*b45.*b13-b13.*b55.*b13.*b44-b11.*b44.*b35.*b35-b14.*b35.*b45.*b13-b45.*b11.*b45.*b33+b15.*b35.*b13.*b44+b55.*b14.*b13.*b34-b55.*b11.*b34.*b34+b45.*b15.*b33.*b14+b34.*b15.*b15.*b34-b34.*b14.*b15.*b35+b11.*b34.*b45.*b35-b45.*b34.*b15.*b13-b44.*b15.*b15.*b33+b11.*b45.*b35.*b34+b55.*b11.*b33.*b44)./d;
ib23=-(b14.*b23.*b45.*b15-b14.*b23.*b55.*b14+b14.*b55.*b34.*b12-b14.*b34.*b25.*b15-b14.*b45.*b35.*b12+b14.*b35.*b25.*b14+b23.*b55.*b11.*b44-b23.*b44.*b15.*b15+b23.*b45.*b15.*b14-b23.*b45.*b11.*b45+b45.*b11.*b35.*b24+b34.*b25.*b11.*b45+b44.*b15.*b35.*b12+b44.*b13.*b25.*b15-b55.*b12.*b13.*b44-b35.*b24.*b15.*b14-b45.*b13.*b25.*b14-b55.*b34.*b24.*b11-b45.*b15.*b13.*b24-b45.*b15.*b34.*b12+b34.*b24.*b15.*b15+b13.*b24.*b55.*b14-b35.*b25.*b11.*b44+b45.*b45.*b12.*b13)./d;
ib24=(-b14.*b23.*b13.*b55+b14.*b23.*b15.*b35+b14.*b13.*b35.*b25+b14.*b33.*b12.*b55-b14.*b35.*b12.*b35-b14.*b15.*b33.*b25+b23.*b45.*b13.*b15+b23.*b11.*b34.*b55-b23.*b11.*b45.*b35-b23.*b34.*b15.*b15+b13.*b13.*b24.*b55-b45.*b13.*b13.*b25-b35.*b24.*b13.*b15+b34.*b15.*b13.*b25+b11.*b35.*b24.*b35+b34.*b35.*b12.*b15+b45.*b12.*b13.*b35-b11.*b33.*b24.*b55+b33.*b24.*b15.*b15+b11.*b45.*b33.*b25-b15.*b13.*b24.*b35-b45.*b33.*b12.*b15-b34.*b12.*b13.*b55-b35.*b25.*b11.*b34)./d;
ib25=-(-b23.*b13.*b14.*b45+b23.*b14.*b14.*b35+b13.*b14.*b34.*b25+b33.*b12.*b14.*b45-b14.*b14.*b33.*b25-b34.*b12.*b14.*b35+b23.*b44.*b13.*b15-b23.*b11.*b44.*b35+b23.*b11.*b34.*b45-b23.*b34.*b14.*b15+b13.*b13.*b24.*b45-b34.*b24.*b13.*b15-b44.*b13.*b13.*b25-b11.*b33.*b24.*b45+b33.*b24.*b14.*b15-b34.*b12.*b13.*b45-b44.*b33.*b12.*b15+b44.*b12.*b13.*b35+b11.*b44.*b33.*b25-b14.*b13.*b24.*b35+b34.*b14.*b13.*b25+b11.*b34.*b24.*b35-b11.*b34.*b34.*b25+b34.*b34.*b12.*b15)./d;
ib33=-(b14.*b24.*b25.*b15-b14.*b24.*b55.*b12+b14.*b55.*b22.*b14-b14.*b22.*b45.*b15+b14.*b25.*b45.*b12-b14.*b25.*b14.*b25-b24.*b24.*b15.*b15+b24.*b45.*b15.*b12+b24.*b55.*b11.*b24-b24.*b25.*b11.*b45-b55.*b11.*b22.*b44-b45.*b12.*b45.*b12+b24.*b12.*b45.*b15+b25.*b45.*b14.*b12+b22.*b45.*b11.*b45-b44.*b25.*b15.*b12+b44.*b12.*b55.*b12+b24.*b15.*b14.*b25+b44.*b22.*b15.*b15-b25.*b44.*b15.*b12-b45.*b15.*b22.*b14+b25.*b11.*b44.*b25-b25.*b45.*b11.*b24-b55.*b14.*b12.*b24)./d;
ib34=-(-b14.*b13.*b55.*b22+b14.*b13.*b25.*b25-b14.*b12.*b35.*b25+b14.*b15.*b35.*b22+b14.*b55.*b23.*b12-b14.*b25.*b23.*b15-b13.*b25.*b12.*b45+b13.*b45.*b15.*b22-b13.*b24.*b15.*b25+b13.*b55.*b12.*b24-b11.*b45.*b35.*b22+b12.*b35.*b12.*b45+b55.*b11.*b22.*b34+b25.*b15.*b12.*b34-b55.*b11.*b23.*b24+b25.*b11.*b23.*b45-b45.*b15.*b23.*b12-b55.*b12.*b12.*b34-b34.*b15.*b15.*b22+b34.*b12.*b15.*b25-b15.*b35.*b12.*b24+b24.*b15.*b23.*b15+b11.*b24.*b35.*b25-b11.*b34.*b25.*b25)./d;
ib35=(b24.*b13.*b14.*b25-b12.*b14.*b24.*b35-b22.*b13.*b14.*b45+b23.*b12.*b14.*b45+b14.*b14.*b22.*b35-b23.*b14.*b14.*b25-b24.*b24.*b13.*b15-b11.*b24.*b34.*b25+b24.*b34.*b12.*b15+b11.*b24.*b24.*b35-b44.*b13.*b12.*b25+b24.*b13.*b12.*b45+b44.*b22.*b13.*b15-b23.*b44.*b12.*b15-b11.*b44.*b22.*b35-b22.*b34.*b14.*b15-b24.*b14.*b12.*b35+b11.*b23.*b44.*b25-b12.*b12.*b34.*b45+b14.*b12.*b34.*b25-b11.*b23.*b24.*b45+b23.*b24.*b14.*b15+b11.*b22.*b34.*b45+b44.*b12.*b12.*b35)./d;
ib44=(-b23.*b13.*b15.*b25+b23.*b13.*b55.*b12-b23.*b15.*b35.*b12+b23.*b11.*b35.*b25+b23.*b15.*b23.*b15-b23.*b55.*b11.*b23-b13.*b55.*b22.*b13-b13.*b35.*b25.*b12+b13.*b15.*b22.*b35+b13.*b13.*b25.*b25-b11.*b33.*b25.*b25-b15.*b22.*b33.*b15+b15.*b33.*b25.*b12+b55.*b11.*b22.*b33+b15.*b33.*b25.*b12+b55.*b23.*b12.*b13-b33.*b12.*b55.*b12+b35.*b12.*b35.*b12+b35.*b25.*b11.*b23-b11.*b35.*b22.*b35-b12.*b13.*b35.*b25+b15.*b35.*b22.*b13-b13.*b25.*b23.*b15-b15.*b35.*b23.*b12)./d;
ib45=(b24.*b23.*b13.*b15-b11.*b24.*b23.*b35-b13.*b12.*b23.*b45-b23.*b23.*b14.*b15+b14.*b12.*b23.*b35+b11.*b23.*b23.*b45-b13.*b13.*b24.*b25-b33.*b24.*b12.*b15+b11.*b33.*b24.*b25+b24.*b12.*b13.*b35-b22.*b34.*b13.*b15+b22.*b13.*b13.*b45+b34.*b13.*b12.*b25-b34.*b12.*b12.*b35+b22.*b33.*b14.*b15+b34.*b23.*b12.*b15-b11.*b22.*b33.*b45+b14.*b13.*b23.*b25+b33.*b12.*b12.*b45-b11.*b34.*b23.*b25+b11.*b22.*b34.*b35-b23.*b12.*b13.*b45-b33.*b14.*b12.*b25-b22.*b14.*b13.*b35)./d;
ib55=-(-b11.*b22.*b33.*b44+b11.*b22.*b34.*b34-b11.*b34.*b23.*b24+b11.*b33.*b24.*b24-b11.*b24.*b23.*b34+b11.*b23.*b23.*b44+b22.*b13.*b13.*b44-b22.*b34.*b13.*b14-b22.*b14.*b13.*b34+b22.*b33.*b14.*b14+b14.*b13.*b23.*b24+b14.*b12.*b23.*b34-b34.*b12.*b12.*b34+b34.*b13.*b12.*b24-b23.*b23.*b14.*b14-b33.*b24.*b12.*b14-b33.*b14.*b12.*b24+b33.*b12.*b12.*b44-b13.*b12.*b23.*b44-b23.*b12.*b13.*b44-b13.*b13.*b24.*b24+b24.*b23.*b13.*b14+b34.*b23.*b12.*b14+b24.*b12.*b13.*b34)./d;
%}
d = (b11.*b22.*b33.*b44.*b55 - b11.*b22.*b33.*b45.*b54 - b11.*b22.*b34.*b43.*b55 + b11.*b22.*b34.*b45.*b53 + b11.*b22.*b35.*b43.*b54 - b11.*b22.*b35.*b44.*b53 - b11.*b23.*b32.*b44.*b55 + b11.*b23.*b32.*b45.*b54 + b11.*b23.*b34.*b42.*b55 - b11.*b23.*b34.*b45.*b52 - b11.*b23.*b35.*b42.*b54 + b11.*b23.*b35.*b44.*b52 + b11.*b24.*b32.*b43.*b55 - b11.*b24.*b32.*b45.*b53 - b11.*b24.*b33.*b42.*b55 + b11.*b24.*b33.*b45.*b52 + b11.*b24.*b35.*b42.*b53 - b11.*b24.*b35.*b43.*b52 - b11.*b25.*b32.*b43.*b54 + b11.*b25.*b32.*b44.*b53 + b11.*b25.*b33.*b42.*b54 - b11.*b25.*b33.*b44.*b52 - b11.*b25.*b34.*b42.*b53 + b11.*b25.*b34.*b43.*b52 - b12.*b21.*b33.*b44.*b55 + b12.*b21.*b33.*b45.*b54 + b12.*b21.*b34.*b43.*b55 - b12.*b21.*b34.*b45.*b53 - b12.*b21.*b35.*b43.*b54 + b12.*b21.*b35.*b44.*b53...
    + b12.*b23.*b31.*b44.*b55 - b12.*b23.*b31.*b45.*b54 - b12.*b23.*b34.*b41.*b55 + b12.*b23.*b34.*b45.*b51 + b12.*b23.*b35.*b41.*b54 - b12.*b23.*b35.*b44.*b51 - b12.*b24.*b31.*b43.*b55 + b12.*b24.*b31.*b45.*b53 + b12.*b24.*b33.*b41.*b55 - b12.*b24.*b33.*b45.*b51 - b12.*b24.*b35.*b41.*b53 + b12.*b24.*b35.*b43.*b51 + b12.*b25.*b31.*b43.*b54 - b12.*b25.*b31.*b44.*b53 - b12.*b25.*b33.*b41.*b54 + b12.*b25.*b33.*b44.*b51 + b12.*b25.*b34.*b41.*b53 - b12.*b25.*b34.*b43.*b51 + b13.*b21.*b32.*b44.*b55 - b13.*b21.*b32.*b45.*b54 - b13.*b21.*b34.*b42.*b55 + b13.*b21.*b34.*b45.*b52 + b13.*b21.*b35.*b42.*b54 - b13.*b21.*b35.*b44.*b52 - b13.*b22.*b31.*b44.*b55 + b13.*b22.*b31.*b45.*b54 + b13.*b22.*b34.*b41.*b55 - b13.*b22.*b34.*b45.*b51 - b13.*b22.*b35.*b41.*b54 + b13.*b22.*b35.*b44.*b51...
    + b13.*b24.*b31.*b42.*b55 - b13.*b24.*b31.*b45.*b52 - b13.*b24.*b32.*b41.*b55 + b13.*b24.*b32.*b45.*b51 + b13.*b24.*b35.*b41.*b52 - b13.*b24.*b35.*b42.*b51 - b13.*b25.*b31.*b42.*b54 + b13.*b25.*b31.*b44.*b52 + b13.*b25.*b32.*b41.*b54 - b13.*b25.*b32.*b44.*b51 - b13.*b25.*b34.*b41.*b52 + b13.*b25.*b34.*b42.*b51 - b14.*b21.*b32.*b43.*b55 + b14.*b21.*b32.*b45.*b53 + b14.*b21.*b33.*b42.*b55 - b14.*b21.*b33.*b45.*b52 - b14.*b21.*b35.*b42.*b53 + b14.*b21.*b35.*b43.*b52 + b14.*b22.*b31.*b43.*b55 - b14.*b22.*b31.*b45.*b53 - b14.*b22.*b33.*b41.*b55 + b14.*b22.*b33.*b45.*b51 + b14.*b22.*b35.*b41.*b53 - b14.*b22.*b35.*b43.*b51 - b14.*b23.*b31.*b42.*b55 + b14.*b23.*b31.*b45.*b52 + b14.*b23.*b32.*b41.*b55 - b14.*b23.*b32.*b45.*b51 - b14.*b23.*b35.*b41.*b52 + b14.*b23.*b35.*b42.*b51 + b14.*b25.*b31.*b42.*b53 - b14.*b25.*b31.*b43.*b52 - b14.*b25.*b32.*b41.*b53 + b14.*b25.*b32.*b43.*b51 + b14.*b25.*b33.*b41.*b52 - b14.*b25.*b33.*b42.*b51 + b15.*b21.*b32.*b43.*b54 - b15.*b21.*b32.*b44.*b53 - b15.*b21.*b33.*b42.*b54 + b15.*b21.*b33.*b44.*b52 + b15.*b21.*b34.*b42.*b53 - b15.*b21.*b34.*b43.*b52 - b15.*b22.*b31.*b43.*b54 + b15.*b22.*b31.*b44.*b53 + b15.*b22.*b33.*b41.*b54 - b15.*b22.*b33.*b44.*b51 - b15.*b22.*b34.*b41.*b53 + b15.*b22.*b34.*b43.*b51 + b15.*b23.*b31.*b42.*b54 - b15.*b23.*b31.*b44.*b52 - b15.*b23.*b32.*b41.*b54 + b15.*b23.*b32.*b44.*b51 + b15.*b23.*b34.*b41.*b52 - b15.*b23.*b34.*b42.*b51 - b15.*b24.*b31.*b42.*b53 + b15.*b24.*b31.*b43.*b52 + b15.*b24.*b32.*b41.*b53 - b15.*b24.*b32.*b43.*b51 - b15.*b24.*b33.*b41.*b52 + b15.*b24.*b33.*b42.*b51);
ib11 = (b22.*b33.*b44.*b55 - b22.*b33.*b45.*b54 - b22.*b34.*b43.*b55 + b22.*b34.*b45.*b53 + b22.*b35.*b43.*b54 - b22.*b35.*b44.*b53 - b23.*b32.*b44.*b55 + b23.*b32.*b45.*b54 + b23.*b34.*b42.*b55 - b23.*b34.*b45.*b52 - b23.*b35.*b42.*b54 + b23.*b35.*b44.*b52 + b24.*b32.*b43.*b55 - b24.*b32.*b45.*b53 - b24.*b33.*b42.*b55 + b24.*b33.*b45.*b52 + b24.*b35.*b42.*b53 - b24.*b35.*b43.*b52 - b25.*b32.*b43.*b54 + b25.*b32.*b44.*b53 + b25.*b33.*b42.*b54 - b25.*b33.*b44.*b52 - b25.*b34.*b42.*b53 + b25.*b34.*b43.*b52)./d;
ib12 = -(b12.*b33.*b44.*b55 - b12.*b33.*b45.*b54 - b12.*b34.*b43.*b55 + b12.*b34.*b45.*b53 + b12.*b35.*b43.*b54 - b12.*b35.*b44.*b53 - b13.*b32.*b44.*b55 + b13.*b32.*b45.*b54 + b13.*b34.*b42.*b55 - b13.*b34.*b45.*b52 - b13.*b35.*b42.*b54 + b13.*b35.*b44.*b52 + b14.*b32.*b43.*b55 - b14.*b32.*b45.*b53 - b14.*b33.*b42.*b55 + b14.*b33.*b45.*b52 + b14.*b35.*b42.*b53 - b14.*b35.*b43.*b52 - b15.*b32.*b43.*b54 + b15.*b32.*b44.*b53 + b15.*b33.*b42.*b54 - b15.*b33.*b44.*b52 - b15.*b34.*b42.*b53 + b15.*b34.*b43.*b52)./d;
ib13 = (b12.*b23.*b44.*b55 - b12.*b23.*b45.*b54 - b12.*b24.*b43.*b55 + b12.*b24.*b45.*b53 + b12.*b25.*b43.*b54 - b12.*b25.*b44.*b53 - b13.*b22.*b44.*b55 + b13.*b22.*b45.*b54 + b13.*b24.*b42.*b55 - b13.*b24.*b45.*b52 - b13.*b25.*b42.*b54 + b13.*b25.*b44.*b52 + b14.*b22.*b43.*b55 - b14.*b22.*b45.*b53 - b14.*b23.*b42.*b55 + b14.*b23.*b45.*b52 + b14.*b25.*b42.*b53 - b14.*b25.*b43.*b52 - b15.*b22.*b43.*b54 + b15.*b22.*b44.*b53 + b15.*b23.*b42.*b54 - b15.*b23.*b44.*b52 - b15.*b24.*b42.*b53 + b15.*b24.*b43.*b52)./d;
ib14 = -(b12.*b23.*b34.*b55 - b12.*b23.*b35.*b54 - b12.*b24.*b33.*b55 + b12.*b24.*b35.*b53 + b12.*b25.*b33.*b54 - b12.*b25.*b34.*b53 - b13.*b22.*b34.*b55 + b13.*b22.*b35.*b54 + b13.*b24.*b32.*b55 - b13.*b24.*b35.*b52 - b13.*b25.*b32.*b54 + b13.*b25.*b34.*b52 + b14.*b22.*b33.*b55 - b14.*b22.*b35.*b53 - b14.*b23.*b32.*b55 + b14.*b23.*b35.*b52 + b14.*b25.*b32.*b53 - b14.*b25.*b33.*b52 - b15.*b22.*b33.*b54 + b15.*b22.*b34.*b53 + b15.*b23.*b32.*b54 - b15.*b23.*b34.*b52 - b15.*b24.*b32.*b53 + b15.*b24.*b33.*b52)./d;
ib15 = (b12.*b23.*b34.*b45 - b12.*b23.*b35.*b44 - b12.*b24.*b33.*b45 + b12.*b24.*b35.*b43 + b12.*b25.*b33.*b44 - b12.*b25.*b34.*b43 - b13.*b22.*b34.*b45 + b13.*b22.*b35.*b44 + b13.*b24.*b32.*b45 - b13.*b24.*b35.*b42 - b13.*b25.*b32.*b44 + b13.*b25.*b34.*b42 + b14.*b22.*b33.*b45 - b14.*b22.*b35.*b43 - b14.*b23.*b32.*b45 + b14.*b23.*b35.*b42 + b14.*b25.*b32.*b43 - b14.*b25.*b33.*b42 - b15.*b22.*b33.*b44 + b15.*b22.*b34.*b43 + b15.*b23.*b32.*b44 - b15.*b23.*b34.*b42 - b15.*b24.*b32.*b43 + b15.*b24.*b33.*b42)./d;
ib22 = (b11.*b33.*b44.*b55 - b11.*b33.*b45.*b54 - b11.*b34.*b43.*b55 + b11.*b34.*b45.*b53 + b11.*b35.*b43.*b54 - b11.*b35.*b44.*b53 - b13.*b31.*b44.*b55 + b13.*b31.*b45.*b54 + b13.*b34.*b41.*b55 - b13.*b34.*b45.*b51 - b13.*b35.*b41.*b54 + b13.*b35.*b44.*b51 + b14.*b31.*b43.*b55 - b14.*b31.*b45.*b53 - b14.*b33.*b41.*b55 + b14.*b33.*b45.*b51 + b14.*b35.*b41.*b53 - b14.*b35.*b43.*b51 - b15.*b31.*b43.*b54 + b15.*b31.*b44.*b53 + b15.*b33.*b41.*b54 - b15.*b33.*b44.*b51 - b15.*b34.*b41.*b53 + b15.*b34.*b43.*b51)./d;
ib23 = -(b11.*b23.*b44.*b55 - b11.*b23.*b45.*b54 - b11.*b24.*b43.*b55 + b11.*b24.*b45.*b53 + b11.*b25.*b43.*b54 - b11.*b25.*b44.*b53 - b13.*b21.*b44.*b55 + b13.*b21.*b45.*b54 + b13.*b24.*b41.*b55 - b13.*b24.*b45.*b51 - b13.*b25.*b41.*b54 + b13.*b25.*b44.*b51 + b14.*b21.*b43.*b55 - b14.*b21.*b45.*b53 - b14.*b23.*b41.*b55 + b14.*b23.*b45.*b51 + b14.*b25.*b41.*b53 - b14.*b25.*b43.*b51 - b15.*b21.*b43.*b54 + b15.*b21.*b44.*b53 + b15.*b23.*b41.*b54 - b15.*b23.*b44.*b51 - b15.*b24.*b41.*b53 + b15.*b24.*b43.*b51)./d;
ib24 = (b11.*b23.*b34.*b55 - b11.*b23.*b35.*b54 - b11.*b24.*b33.*b55 + b11.*b24.*b35.*b53 + b11.*b25.*b33.*b54 - b11.*b25.*b34.*b53 - b13.*b21.*b34.*b55 + b13.*b21.*b35.*b54 + b13.*b24.*b31.*b55 - b13.*b24.*b35.*b51 - b13.*b25.*b31.*b54 + b13.*b25.*b34.*b51 + b14.*b21.*b33.*b55 - b14.*b21.*b35.*b53 - b14.*b23.*b31.*b55 + b14.*b23.*b35.*b51 + b14.*b25.*b31.*b53 - b14.*b25.*b33.*b51 - b15.*b21.*b33.*b54 + b15.*b21.*b34.*b53 + b15.*b23.*b31.*b54 - b15.*b23.*b34.*b51 - b15.*b24.*b31.*b53 + b15.*b24.*b33.*b51)./d;
ib25 = -(b11.*b23.*b34.*b45 - b11.*b23.*b35.*b44 - b11.*b24.*b33.*b45 + b11.*b24.*b35.*b43 + b11.*b25.*b33.*b44 - b11.*b25.*b34.*b43 - b13.*b21.*b34.*b45 + b13.*b21.*b35.*b44 + b13.*b24.*b31.*b45 - b13.*b24.*b35.*b41 - b13.*b25.*b31.*b44 + b13.*b25.*b34.*b41 + b14.*b21.*b33.*b45 - b14.*b21.*b35.*b43 - b14.*b23.*b31.*b45 + b14.*b23.*b35.*b41 + b14.*b25.*b31.*b43 - b14.*b25.*b33.*b41 - b15.*b21.*b33.*b44 + b15.*b21.*b34.*b43 + b15.*b23.*b31.*b44 - b15.*b23.*b34.*b41 - b15.*b24.*b31.*b43 + b15.*b24.*b33.*b41)./d;
ib33 = (b11.*b22.*b44.*b55 - b11.*b22.*b45.*b54 - b11.*b24.*b42.*b55 + b11.*b24.*b45.*b52 + b11.*b25.*b42.*b54 - b11.*b25.*b44.*b52 - b12.*b21.*b44.*b55 + b12.*b21.*b45.*b54 + b12.*b24.*b41.*b55 - b12.*b24.*b45.*b51 - b12.*b25.*b41.*b54 + b12.*b25.*b44.*b51 + b14.*b21.*b42.*b55 - b14.*b21.*b45.*b52 - b14.*b22.*b41.*b55 + b14.*b22.*b45.*b51 + b14.*b25.*b41.*b52 - b14.*b25.*b42.*b51 - b15.*b21.*b42.*b54 + b15.*b21.*b44.*b52 + b15.*b22.*b41.*b54 - b15.*b22.*b44.*b51 - b15.*b24.*b41.*b52 + b15.*b24.*b42.*b51)./d;
ib34 = -(b11.*b22.*b34.*b55 - b11.*b22.*b35.*b54 - b11.*b24.*b32.*b55 + b11.*b24.*b35.*b52 + b11.*b25.*b32.*b54 - b11.*b25.*b34.*b52 - b12.*b21.*b34.*b55 + b12.*b21.*b35.*b54 + b12.*b24.*b31.*b55 - b12.*b24.*b35.*b51 - b12.*b25.*b31.*b54 + b12.*b25.*b34.*b51 + b14.*b21.*b32.*b55 - b14.*b21.*b35.*b52 - b14.*b22.*b31.*b55 + b14.*b22.*b35.*b51 + b14.*b25.*b31.*b52 - b14.*b25.*b32.*b51 - b15.*b21.*b32.*b54 + b15.*b21.*b34.*b52 + b15.*b22.*b31.*b54 - b15.*b22.*b34.*b51 - b15.*b24.*b31.*b52 + b15.*b24.*b32.*b51)./d;
ib35 = (b11.*b22.*b34.*b45 - b11.*b22.*b35.*b44 - b11.*b24.*b32.*b45 + b11.*b24.*b35.*b42 + b11.*b25.*b32.*b44 - b11.*b25.*b34.*b42 - b12.*b21.*b34.*b45 + b12.*b21.*b35.*b44 + b12.*b24.*b31.*b45 - b12.*b24.*b35.*b41 - b12.*b25.*b31.*b44 + b12.*b25.*b34.*b41 + b14.*b21.*b32.*b45 - b14.*b21.*b35.*b42 - b14.*b22.*b31.*b45 + b14.*b22.*b35.*b41 + b14.*b25.*b31.*b42 - b14.*b25.*b32.*b41 - b15.*b21.*b32.*b44 + b15.*b21.*b34.*b42 + b15.*b22.*b31.*b44 - b15.*b22.*b34.*b41 - b15.*b24.*b31.*b42 + b15.*b24.*b32.*b41)./d;
ib44 = (b11.*b22.*b33.*b55 - b11.*b22.*b35.*b53 - b11.*b23.*b32.*b55 + b11.*b23.*b35.*b52 + b11.*b25.*b32.*b53 - b11.*b25.*b33.*b52 - b12.*b21.*b33.*b55 + b12.*b21.*b35.*b53 + b12.*b23.*b31.*b55 - b12.*b23.*b35.*b51 - b12.*b25.*b31.*b53 + b12.*b25.*b33.*b51 + b13.*b21.*b32.*b55 - b13.*b21.*b35.*b52 - b13.*b22.*b31.*b55 + b13.*b22.*b35.*b51 + b13.*b25.*b31.*b52 - b13.*b25.*b32.*b51 - b15.*b21.*b32.*b53 + b15.*b21.*b33.*b52 + b15.*b22.*b31.*b53 - b15.*b22.*b33.*b51 - b15.*b23.*b31.*b52 + b15.*b23.*b32.*b51)./d;
ib45 = -(b11.*b22.*b33.*b45 - b11.*b22.*b35.*b43 - b11.*b23.*b32.*b45 + b11.*b23.*b35.*b42 + b11.*b25.*b32.*b43 - b11.*b25.*b33.*b42 - b12.*b21.*b33.*b45 + b12.*b21.*b35.*b43 + b12.*b23.*b31.*b45 - b12.*b23.*b35.*b41 - b12.*b25.*b31.*b43 + b12.*b25.*b33.*b41 + b13.*b21.*b32.*b45 - b13.*b21.*b35.*b42 - b13.*b22.*b31.*b45 + b13.*b22.*b35.*b41 + b13.*b25.*b31.*b42 - b13.*b25.*b32.*b41 - b15.*b21.*b32.*b43 + b15.*b21.*b33.*b42 + b15.*b22.*b31.*b43 - b15.*b22.*b33.*b41 - b15.*b23.*b31.*b42 + b15.*b23.*b32.*b41)./d;
ib55 = (b11.*b22.*b33.*b44 - b11.*b22.*b34.*b43 - b11.*b23.*b32.*b44 + b11.*b23.*b34.*b42 + b11.*b24.*b32.*b43 - b11.*b24.*b33.*b42 - b12.*b21.*b33.*b44 + b12.*b21.*b34.*b43 + b12.*b23.*b31.*b44 - b12.*b23.*b34.*b41 - b12.*b24.*b31.*b43 + b12.*b24.*b33.*b41 + b13.*b21.*b32.*b44 - b13.*b21.*b34.*b42 - b13.*b22.*b31.*b44 + b13.*b22.*b34.*b41 + b13.*b24.*b31.*b42 - b13.*b24.*b32.*b41 - b14.*b21.*b32.*b43 + b14.*b21.*b33.*b42 + b14.*b22.*b31.*b43 - b14.*b22.*b33.*b41 - b14.*b23.*b31.*b42 + b14.*b23.*b32.*b41)./d;
% y project onto B'
py1 = sum(conj(col1).*y,1);
py2 = sum(conj(col2).*y,1);
py3 = sum(conj(col3).*y,1);
py4 = sum(conj(col4).*y,1);
py5 = sum(conj(col5).*y,1);
% calculate x
x(1,:) = sum(ib11.*py1 + ib12.*py2 + ib13.*py3 + ib14.*py4 + ib15.*py5,1);
x(2,:) = sum(conj(ib12).*py1 + ib22.*py2 + ib23.*py3 + ib24.*py4 + ib25.*py5,1);
x(3,:) = sum(conj(ib13).*py1 + conj(ib23).*py2 + ib33.*py3 + ib34.*py4 + ib35.*py5,1);
x(4,:) = sum(conj(ib14).*py1 + conj(ib24).*py2 + conj(ib34).*py3 + ib44.*py4 + ib45.*py5,1);
x(5,:) = sum(conj(ib15).*py1 + conj(ib25).*py2 + conj(ib35).*py3 + conj(ib45).*py4 + ib55.*py5,1);
end
%{
function x = matrixB(col1, col2, col3, col4, col5)
% assemble B'*B
b11 = sum(conj(col1).*col1,1);x(1,1,:) = b11;
b12 = sum(conj(col1).*col2,1);x(1,2,:) = b12;x(2,1,:) = conj(b12);
b13 = sum(conj(col1).*col3,1);x(1,3,:) = b13;x(3,1,:) = conj(b13);
b14 = sum(conj(col1).*col4,1);x(1,4,:) = b14;x(4,1,:) = conj(b14);
b15 = sum(conj(col1).*col5,1);x(1,5,:) = b15;x(5,1,:) = conj(b15);
b22 = sum(conj(col2).*col2,1);x(2,2,:) = b22;
b23 = sum(conj(col2).*col3,1);x(2,3,:) = b23;x(3,2,:) = conj(b23);
b24 = sum(conj(col2).*col4,1);x(2,4,:) = b24;x(4,2,:) = conj(b24);
b25 = sum(conj(col2).*col5,1);x(2,5,:) = b25;x(5,2,:) = conj(b25);
b33 = sum(conj(col3).*col3,1);x(3,3,:) = b33;
b34 = sum(conj(col3).*col4,1);x(3,4,:) = b34;x(4,3,:) = conj(b34);
b35 = sum(conj(col3).*col5,1);x(3,5,:) = b35;x(5,3,:) = conj(b35);
b44 = sum(conj(col4).*col4,1);x(4,4,:) = b44;
b45 = sum(conj(col4).*col5,1);x(4,5,:) = b45;x(5,4,:) = conj(b45);
b55 = sum(conj(col5).*col5,1);x(5,5,:) = b55;
% x = [b11,b12,b13,b14,b15; ...
%     conj(b12),b22,b23,b24,b25; ...
%     conj(b13),conj(b23),b33,b34,b35;...
%     conj(b14),conj(b24),conj(b34),b44,b35;...
%     conj(b15),conj(b25),conj(b35),b45,b55];
end
%}