% clear all;
% x:horizontal axis, y:vertical axis,

tic
clear
close all
DataPath = 'D:\WorkingData\FaceRecog_VP';

% l = [253, 132, 259, 210; 373 18 371 126; 541 50 522 197];
% img = double(imread(fullfile(DataPath,'.\2_640x480\C1_000000328.bmp')))/255;
% l = [320, 373, 338, 481; 34 560 80 665; 1149 279 1079 498];
% img = double(imread('pic.JPG'))/255;
% l = [111, 95, 172, 412; 383 150 379 274; 482 189 467 306];
% img = double(imread(fullfile(DataPath,'Picture 15.JPG')))/255;
% l = [227, 133, 251, 358; 513 134 497 346];
% img = double(imread('sync01.png'))/255;
% l = [80, 197, 122, 368; 207 212 220 312; 598 137 550 330];
% img = double(imread(fullfile(DataPath,'sync02.png')))/255;
% l = [26, 16, 58, 136; 3 38 38 165; 155 4 165 56];
% img = double(imread(fullfile(DataPath,'C1_000012440.bmp')))/255;
[img, l] = genPerspective01();
h = size(img,1);
w = size(img,2);

% calculate VP
for i = 1:size(l,1);
    % x = ay+b;
    a(i) = (l(i, 3) - l(i,1)) / (l(i,4) - l(i,2));
    b(i) = l(i,1) - l(i,2) * a(i);
end

% x - y*a1 = b1 ...
% [1 -a1] [x, y]' = b1
VP = [ ones(length(a),1) -a'] \ b';     % top-left

IC = [(w+1)/2; (h+1)/2];                % top-left

img_line = img;
for x = 10:30:size(img_line,2)-1
    a_tmp = (x - VP(1)) / (1 - VP(2));
    b_tmp = x - a_tmp;
    
    for j = 1:size(img_line,1);
        img_line( j, round(a_tmp*j + b_tmp), 1) = 255;
        img_line( j, round(a_tmp*j + b_tmp), 2) = 0;
        img_line( j, round(a_tmp*j + b_tmp), 3) = 0;
    end
end


% parameter of central separation line
dPVIC = VP-IC;
dPVIC_normal = normc(dPVIC);

% mRotateInv * VP and mRotateInv * IC should on the same vertical line
mRotateInv = [dPVIC_normal(2) -dPVIC_normal(1);...
    dPVIC_normal(1) dPVIC_normal(2)];
mRotate = [dPVIC_normal(2) dPVIC_normal(1);...
    -dPVIC_normal(1) dPVIC_normal(2)];

% find size of R_I image
tx = repmat(1:w, h, 1);
tx = tx(:)';
ty = repmat(1:h, 1, w);
img_pos_c = [tx - IC(1); ty - IC(2)];

tmp = mRotateInv*img_pos_c;
x_m = min(tmp(1,:));
x_M = max(tmp(1,:));
y_m = min(tmp(2,:));
y_M = max(tmp(2,:));

img_R_I_xy2ijShift = ceil([max(abs(tmp(1,:))); ...
    max(abs(tmp(2,:)))]) + [1;1];

img_R_I_w = img_R_I_xy2ijShift(1) * 2 - 1;
img_R_I_h = img_R_I_xy2ijShift(2) * 2 - 1;

img_R_I_IC = [0;0] + img_R_I_xy2ijShift;
img_R_I_VP = mRotateInv*(VP-IC) + img_R_I_xy2ijShift;


% <-- rectify rotation
tx = repmat( ((1:img_R_I_w) - img_R_I_xy2ijShift(1)), img_R_I_h, 1);
tx = tx(:)';
ty = repmat( ((1:img_R_I_h) - img_R_I_xy2ijShift(2)), 1, img_R_I_w);

tmp = mRotate*[tx; ty];
tmp(1,:) = tmp(1,:) + IC(1);
tmp(2,:) = tmp(2,:) + IC(2);

img_R_I = plotImgPoint(img, tmp, [img_R_I_w, img_R_I_h]);
img_R_I_line = plotImgPoint(img_line, tmp, [img_R_I_w, img_R_I_h]);

% --> <End> rectify rotation





tmpx = zeros(img_R_I_h, img_R_I_w);
tmpy = zeros(img_R_I_h, img_R_I_w);

% now suppose the center is on img_rotateInv's top-middle
zeroCenter = [0; 0];
shiftToZeroCenter = [zeroCenter(1)-img_R_I_IC(1); 0];
IC_zeroCenter = img_R_I_IC + shiftToZeroCenter;
VP_zeroCenter = img_R_I_VP + shiftToZeroCenter;

% in fact, I only care about the vertical distance from IC to VP now
% (rotation rectified)



% assum Z and f, from f and VP, we can get theta
Z = 500;
f = 500;

% img_R_I_IC_xy = img_R_I_IC - img_R_I_IC;
% img_R_I_VP_xy = img_R_I_VP - img_R_I_IC;
img_R_I_IC_xy = [0; 0];
img_R_I_VP_xy = mRotateInv * (VP-IC);

VP_dY = img_R_I_VP_xy(2);
theta = atan(f/VP_dY);
theta_ground = theta - pi/2;



img_rect_w = round(img_R_I_w*1.5);
img_rect_h = round(img_R_I_h*5);
img_rect_ij2xyShift = -[(img_rect_w+1)/2; (img_rect_h+1)/2];

% img_rect = zeros(img_rect_h, img_rect_w,3);
img_rect_R_P_idx_X = zeros(img_rect_h, img_rect_w);
img_rect_R_P_idx_Y = zeros(img_rect_h, img_rect_w);
for i = 1: img_rect_w
    x0 = i+img_rect_ij2xyShift(1);
    for j = 1: img_rect_h
        y0 = j+img_rect_ij2xyShift(2);
        
        x1 = x0;
        y1 = cos(theta)*y0;
        z1 = Z + sin(theta)*y0;
        
        x1_ = x1 / z1 * f;
        y1_ = y1 / z1 * f;
        
        img_rect_R_P_idx_X(j,i) = x1_;
        img_rect_R_P_idx_Y(j,i) = y1_;
    end
end
tmp = mRotate * [img_rect_R_P_idx_X(:)';img_rect_R_P_idx_Y(:)'];
tmp(1,:) = tmp(1,:) + IC(1);
tmp(2,:) = tmp(2,:) + IC(2);
img_rect = plotImgPoint(img, tmp, [img_rect_w, img_rect_h]);
img_rect_line = plotImgPoint(img_line, tmp, [img_rect_w, img_rect_h]);




% img_rect = zeros(img_rect_h, img_rect_w,3);
img_rect2_h = img_rect_h;
img_rect2_w = img_rect_w;
img_rect2_R_P_idx_X = zeros(img_rect2_h, img_rect2_w);
img_rect2_R_P_idx_Y = zeros(img_rect2_h, img_rect2_w);


    
for i = 1:img_rect2_w
    x0 = i+img_rect_ij2xyShift(1);
    
    th_tmp = x0 / VP_dY;
    mRot_tmp = [cos(th_tmp) -sin(th_tmp); sin(th_tmp) cos(th_tmp)];
    for j = 1:img_rect2_h
        y0 = j+img_rect_ij2xyShift(2);
        y1 = cos(theta)*y0;
        z1 = Z + sin(theta)*y0;
        
        y1_ = y1 / z1 * f;
        % this is only work for the centeral vertical line
        
        tmp = mRot_tmp * [0; y1_-VP_dY];
        img_rect2_R_P_idx_X(j,i) = tmp(1);
        img_rect2_R_P_idx_Y(j,i) = tmp(2)+VP_dY;
    end
end

tmp = mRotate * [img_rect2_R_P_idx_X(:)';img_rect2_R_P_idx_Y(:)'];
tmp(1,:) = tmp(1,:) + IC(1);
tmp(2,:) = tmp(2,:) + IC(2);
img_rect2 = plotImgPoint(img, tmp, [img_rect2_w, img_rect2_h]);
img_rect2_line = plotImgPoint(img_line, tmp, [img_rect2_w, img_rect2_h]);

imshow(img_R_I);
figure, imshow(img_R_I_line);
figure, imshow(img_rect)
figure, imshow(img_rect_line);
figure, imshow(img_rect2);
figure, imshow(img_rect2_line);








return;




% an older intersting approach
Z = 200;
f = 200;
for i = 1:img_R_I_h

    x1 = 100;
    x2 = x1 * (VP_zeroCenter(2) - i)/VP_zeroCenter(2);
    
    X_bar = x1 * Z / f;
    
    l = ( (i*X_bar/x2)^2 + f^2 * ( X_bar/x2 - X_bar/x1 )^2  )^0.5;

    for j = 1:img_R_I_w
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
for i=1:img_R_I_h
    for j = 1:img_R_I_w
        img__(tmpy(i,j),tmpx(i,j)-minX+1,1:3) = img_rotateInvLine(...
            i, j,:);
    end
end


















