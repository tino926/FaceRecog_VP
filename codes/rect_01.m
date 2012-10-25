% clear all;
tic
% l = [253, 132, 259, 210; 373 18 371 126; 541 50 522 197];
% img = double(imread('.\2_640x480\C1_000000328.bmp'))/255;
% l = [320, 373, 338, 481; 34 560 80 665; 1149 279 1079 498];
% img = double(imread('pic.JPG'))/255;
l = [111, 95, 172, 412; 383 150 379 274; 482 189 467 306];
img = double(imread('Picture 15.JPG'))/255;
% l = [227, 133, 251, 358; 513 134 497 346];
% img = double(imread('syn01.png'))/255;
h = size(img,1);
w = size(img,2);

for i = 1:size(l,1);
    % x = ay+b;
    a(i) = (l(i, 3) - l(i,1)) / (l(i,4) - l(i,2));
    b(i) = l(i,1) - l(i,2) * a(i);
end

% x - y*a1 = b1 ...
% [1 -a1] [x, y]' = b1
X = [ ones(length(a),1) -a'] \ b';

C = [(w+1)/2; (h+1)/2];

% parameter of central separation line
tmp = (C+X)/2;
if abs(C(1)-X(1))>0.01
    % to straignt
    theta = atan( (C(1)-X(1)) / (C(2)-X(2)) );
    mRotate = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    mRotateInv = [cos(-theta) -sin(-theta); sin(-theta) cos(-theta)];
else
    theta = 0;
    y_center = tmp(2);
    mRotate = [1. 0.; 0. 1.];
    mRotateInv = [1. 0.; 0. 1.];
end
X_r = mRotate*X;
C_r = mRotate*C;
tmp = mRotate*tmp;
x_c_r = tmp(1);
y_c_r = tmp(2);

% r: rotate, s: shift
C_r_s = C_r - X_r;
x_c_r_s = x_c_r-X_r(1);
y_c_r_s = y_c_r-X_r(2);


x_rect_min = Inf;
x_rect_max = -Inf;
y_rect_min = Inf;
y_rect_max = -Inf;


xy_b = [ 1:w ones(1,h) 1:w ones(1,h)*w; ...
        ones(1,w) 1:h ones(1,w)*h 1:h];
for i = 1:(h+w)*2
    x_tmp = mRotate*xy_b(:,i);
    x_tmp = x_tmp-X_r;
    
    % x_tmp(2) should not be zero...
    x_rect = x_tmp(1) * y_c_r_s / x_tmp(2);
    y_rect = C_r_s(2) + (x_tmp(2)-C_r_s(2)) * y_c_r_s / x_tmp(2);
    x_rect_min = min(x_rect_min, x_rect);
    x_rect_max = max(x_rect_max, x_rect);
    y_rect_min = min(y_rect_min, y_rect);
    y_rect_max = max(y_rect_max, y_rect);
%     img_rect( ceil(y_rect-y_rect_min), ceil(x_rect-x_rect_min),:) = ...
%             255;
end
x_rect_min = floor(x_rect_min);
x_rect_max = ceil(x_rect_max);
y_rect_min = floor(y_rect_min);
y_rect_max = ceil(y_rect_max);

w_rect = x_rect_max - x_rect_min + 1;
h_rect = y_rect_max - y_rect_min + 1;
img_rect = zeros(h_rect, w_rect,3);

img_rect_direct = img_rect;

for x = 1:w
    for y = 1:h
        x_tmp = mRotate*[x;y];
        x_tmp = x_tmp-X_r;
        
        % x_tmp(2) should not be zero...
        x_rect = x_tmp(1) * y_c_r_s / x_tmp(2);
        y_rect = C_r_s(2) + (x_tmp(2)-C_r_s(2)) * y_c_r_s / x_tmp(2);
        
        img_rect_direct( ceil(y_rect-y_rect_min), ceil(x_rect-x_rect_min),:) = ...
            img(y,x,:);
    end
end


for x = 1:w_rect
    for y = 1:h_rect
        x_t = x + x_rect_min;
        y_t = y + y_rect_min;
        % now x_t = x_rect in above
        y_t = (C_r_s(2) * y_c_r_s) / (C_r_s(2) + y_c_r_s - y_t);
        x_t = x_t*y_t/y_c_r_s;
        x_tmp = [x_t; y_t] + X_r;
        x_tmp = mRotateInv*x_tmp;
        
        x_o_p = floor(x_tmp(1));
        x_o_n = ceil(x_tmp(1));
        y_o_p = floor(x_tmp(2));
        y_o_n = ceil(x_tmp(2));
        
        if (x_o_p > 0 && x_o_n <=w && y_o_p >0 && y_o_n <=h)
%             img_rect(y,x,:) = img(y_o,x_o,:);
            dx_o_p = x_tmp(1) - x_o_p;
            dx_o_n = x_o_n - x_tmp(1);
            dy_o_p = x_tmp(2) - y_o_p;
            dy_o_n = y_o_n - x_tmp(2);
            
            img_rect(y,x,:) = ...
                img(y_o_p, x_o_p, :) * dy_o_n * dx_o_n+ ...
                img(y_o_p, x_o_n, :) * dy_o_n * dx_o_p + ...
                img(y_o_n, x_o_p, :) * dy_o_p * dx_o_n + ...
                img(y_o_n, x_o_n, :) * dy_o_p * dx_o_p;
        end
    end
end












for i = 10:30:size(img,2)-1
    a_tmp = (i - X(1)) / (1 - X(2));
    b_tmp = i - a_tmp;
    
    for j = 1:size(img,1);
        img( j, round(a_tmp*j + b_tmp), 1) = 255;
    end
end

