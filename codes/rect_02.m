% clear all;
% x:horizontal axis, y:vertical axis,

tic
% l = [253, 132, 259, 210; 373 18 371 126; 541 50 522 197];
% img = double(imread('.\2_640x480\C1_000000328.bmp'))/255;
% l = [320, 373, 338, 481; 34 560 80 665; 1149 279 1079 498];
% img = double(imread('pic.JPG'))/255;
% l = [111, 95, 172, 412; 383 150 379 274; 482 189 467 306];
% img = double(imread('Picture 15.JPG'))/255;
% l = [227, 133, 251, 358; 513 134 497 346];
% img = double(imread('sync01.png'))/255;
% l = [80, 197, 122, 368; 207 212 220 312; 598 137 550 330];
% img = double(imread('sync02.png'))/255;
[img, l] = genPerspective01();
h = size(img,1);
w = size(img,2);

for i = 1:size(l,1);
    % x = ay+b;
    a(i) = (l(i, 3) - l(i,1)) / (l(i,4) - l(i,2));
    b(i) = l(i,1) - l(i,2) * a(i);
end

% x - y*a1 = b1 ...
% [1 -a1] [x, y]' = b1
VP = [ ones(length(a),1) -a'] \ b';

IC = [(w+1)/2; (h+1)/2];

% parameter of central separation line
tmp = (IC+VP)/2;

dPVIC = VP-IC;
dPVIC_normal = normc(dPVIC);

% mRotateInv * VP and mRotateInv * IC should on the same vertical line
mRotateInv = [dPVIC_normal(2) -dPVIC_normal(1);...
    dPVIC_normal(1) dPVIC_normal(2)];
mRotate = [dPVIC_normal(2) dPVIC_normal(1);...
    -dPVIC_normal(1) dPVIC_normal(2)];

tmp = [1 w 1 w; 1 1 h h];
tmp = mRotateInv*tmp;

x_m = min(tmp(1,:));
x_M = max(tmp(1,:));
y_m = min(tmp(2,:));
y_M = max(tmp(2,:));

tmp = mRotateInv*IC;

img_rotateInv_shift = zeros(2,1);
img_rotateInv_shift(1) = ceil(tmp(1)+(-x_m)) + 0.5 - tmp(1) + 1;
img_rotateInv_shift(2) = ceil(tmp(2)+(-y_m)) + 0.5 - tmp(2) + 1;

img_rotateInv_w = ceil(x_M + img_rotateInv_shift(1)) + 1;
img_rotateInv_h = ceil(y_M + img_rotateInv_shift(2)) + 1;
img_rotateInv = zeros(img_rotateInv_h, img_rotateInv_w, 3);



IC_rotateInv = mRotateInv*IC + img_rotateInv_shift;
VP_rotateInv = mRotateInv*VP + img_rotateInv_shift;



% <-- rectify rotation
x = repmat( ((1:img_rotateInv_w) - img_rotateInv_shift(1)), ...
    img_rotateInv_h, 1);
x = x(:)';
y = repmat( ((1:img_rotateInv_h) - img_rotateInv_shift(2)), ...
    1, img_rotateInv_w);

tmp = mRotate*[x; y];


x_b = min(max(floor(tmp(1,:)),1),w);
x_u = min(max(ceil(tmp(1,:)),1),w);
x_w_b = 1 - min(abs(tmp(1,:) - x_b), 1);
x_w_u = 1 - x_w_b;

y_b = min(max(floor(tmp(2,:)),1),h);
y_u = min(max(ceil(tmp(2,:)),1),h);
y_w_b = 1 - min(abs(tmp(2,:) - y_b), 1);
y_w_u = 1 - y_w_b;


tmp = floor(((1:length(x_b)*3)-1)/length(x_b)) * w*h;
idx_ybxb = repmat( (x_b-1)*h+y_b , 1, 3) + tmp;
idx_ybxu = repmat( (x_u-1)*h+y_b , 1, 3) + tmp;
idx_yuxb = repmat( (x_b-1)*h+y_u , 1, 3) + tmp;
idx_yuxu = repmat( (x_u-1)*h+y_u , 1, 3) + tmp;
w_ybxb = repmat( x_w_b.*y_w_b , 1, 3);
w_ybxu = repmat( x_w_u.*y_w_b , 1, 3);
w_yuxb = repmat( x_w_b.*y_w_u , 1, 3);
w_yuxu = repmat( x_w_u.*y_w_u , 1, 3);

img_rotateInv(:) = img(idx_ybxb) .* w_ybxb + ...
    img(idx_ybxu) .* w_ybxu + ...
    img(idx_yuxb) .* w_yuxb + ...
    img(idx_yuxu) .* w_yuxu;
imshow(img_rotateInv);
% -->





% img_rotateInvLine = img_rotateInv;
% 
% for i = 10:30:size(img_rotateInvLine,2)-1
%     a_tmp = (i - VP_rotateInv(1)) / (1 - VP_rotateInv(2));
%     b_tmp = i - a_tmp;
%     
%     for j = 1:size(img_rotateInvLine,1);
%         img_rotateInvLine( j, round(a_tmp*j + b_tmp), 1) = 255;
%         img_rotateInvLine( j, round(a_tmp*j + b_tmp), 2) = 0;
%         img_rotateInvLine( j, round(a_tmp*j + b_tmp), 3) = 0;
%     end
% end




tmpx = zeros(img_rotateInv_h, img_rotateInv_w);
tmpy = zeros(img_rotateInv_h, img_rotateInv_w);

% now suppose the center is on img_rotateInv's top-middle
zeroCenter = [0; 0];
shiftToZeroCenter = [zeroCenter(1)-IC_rotateInv(1); 0];
IC_zeroCenter = IC_rotateInv + shiftToZeroCenter;
VP_zeroCenter = VP_rotateInv + shiftToZeroCenter;

img_rotateInvLine = img_rotateInv;

for i = 10:30:size(img_rotateInvLine,2)-1
    x = i + shiftToZeroCenter(1);
    a_tmp = (x - VP_zeroCenter(1)) / (1 - VP_zeroCenter(2));
    b_tmp = x - a_tmp;
    
    for j = 1:size(img_rotateInvLine,1);
        img_rotateInvLine( j, round(a_tmp*j + b_tmp - shiftToZeroCenter(1)), 1) = 255;
        img_rotateInvLine( j, round(a_tmp*j + b_tmp - shiftToZeroCenter(1)), 2) = 0;
        img_rotateInvLine( j, round(a_tmp*j + b_tmp - shiftToZeroCenter(1)), 3) = 0;
%         imshow(img_rotateInvLine);
%         pause(0.1);
    end
end










Z = 10;
f = 20;
for i = 1:img_rotateInv_h

    x1 = 100;
    x2 = x1 * (VP_zeroCenter(2) - i)/VP_zeroCenter(2);
    
    X_bar = x1 * Z / f;
    
    l = ( (i*X_bar/x2)^2 + f^2 * ( X_bar/x2 - X_bar/x1 )^2  )^0.5;

    for j = 1:img_rotateInv_w
        x2 = j+shiftToZeroCenter(1);
        x1 = x2 * VP_zeroCenter(2) / (VP_zeroCenter(2) - i);
        X_bar = x1 * Z / f;
        
        tmpx(i,j) = round(X_bar-shiftToZeroCenter(1));
        tmpy(i,j) = round(l);
    end
end

minX = min(tmpx(:));

img__ = zeros(max(tmpy(:))-min(tmpy(:)) + 1, ...
    max(tmpx(:))-minX+1, 3);
for i=1:img_rotateInv_h
    for j = 1:img_rotateInv_w
        img__(tmpy(i,j),tmpx(i,j)-minX+1,1:3) = img_rotateInvLine(...
            i, j,:);
    end
end






















center_shift = -IC_rotateInv;
IC_center = IC_rotateInv + center_shift;
VP_center = VP_rotateInv + center_shift;

x1 = 10.;
y1 = -1.;

x2 = 20.;
y2 = (y1-VP_center(2)) / (x1-VP_center(1)) * (x2-x1) + y1;



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

