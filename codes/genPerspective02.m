function [p, l] = genPerspective02(params)
% want to generate c-board

T = genTemplate();
T = imresize(T, 1);
[T_H, T_W, nCh] = size(T);


Z = 300;
f = 100;
theta = pi/180 * 30;

% set the center of T to be (0,0)
Tij2xyShift = -[(T_W+1)/2; (T_H+1)/2];

% the pixel on T to where
T_R_P_idx_X = zeros(T_H, T_W); % T->rotate->project   position x
T_R_P_idx_Y = zeros(T_H, T_W);

mRotate = [cos(theta) -sin(theta); sin(theta) cos(theta)];
for i = 1: T_W
    % on T: (x0, y0, Z)
    % on rotated: (x1, y1, z1)
    % on T->R->P: (x1_, y1_, f)
    x0_ = i+Tij2xyShift(1);
    x0 = sin(x0_/Z)*Z;
    z0 = cos(x0_/Z)*Z;
    for j = 1: T_H
        y0 = j+Tij2xyShift(2);
        
        x1 = x0;
%         y1 = cos(theta)*y0;
%         z1 = Z + sin(theta)*y0;
        tmp = mRotate * [y0; z0-Z];
        y1 = tmp(1);
        z1 = tmp(2)+Z;
        
        x1_ = x1 / z1 * f;
        y1_ = y1 / z1 * f;
        
        T_R_P_idx_X(j,i) = x1_;
        T_R_P_idx_Y(j,i) = y1_;
    end
end

% note: 
T_R_P_xy2ijShift = ceil([max(abs(T_R_P_idx_X(:))); ...
    max(abs(T_R_P_idx_Y(:)))]) + [1;1];

img_T_R_P_w = T_R_P_xy2ijShift(1)*2 - 1;
img_T_R_P_h = T_R_P_xy2ijShift(2)*2 - 1;

T_R_P_idx_Xij = round(T_R_P_idx_X + T_R_P_xy2ijShift(1));
T_R_P_idx_Yij = round(T_R_P_idx_Y + T_R_P_xy2ijShift(2));

l = [T_R_P_idx_X(round(T_H*0.25), round(T_W * 0.25)) + T_R_P_xy2ijShift(1), ...
    T_R_P_idx_Y(round(T_H*0.25), round(T_W * 0.25)) + T_R_P_xy2ijShift(2), ...
    T_R_P_idx_X(round(T_H*0.75), round(T_W * 0.25)) + T_R_P_xy2ijShift(1), ...
    T_R_P_idx_Y(round(T_H*0.75), round(T_W * 0.25)) + T_R_P_xy2ijShift(2); ...
    T_R_P_idx_X(round(T_H*0.25), round(T_W * 0.75)) + T_R_P_xy2ijShift(1), ...
    T_R_P_idx_Y(round(T_H*0.25), round(T_W * 0.75)) + T_R_P_xy2ijShift(2), ...
    T_R_P_idx_X(round(T_H*0.75), round(T_W * 0.75)) + T_R_P_xy2ijShift(1), ...
    T_R_P_idx_Y(round(T_H*0.75), round(T_W * 0.75)) + T_R_P_xy2ijShift(2)];

img_T_R_P_1 = zeros(img_T_R_P_h,...
    img_T_R_P_w, nCh);

for i = 1: T_W
    for j = 1: T_H
        img_T_R_P_1( T_R_P_idx_Yij(j,i),...
            T_R_P_idx_Xij(j,i), :) = ...
            T(j,i,:);
    end
end
p = img_T_R_P_1;
return;

% reverse ... 
img_T_R_P_toT_idx_X = zeros(img_T_R_P_h, img_T_R_P_w);
img_T_R_P_toT_idx_Y = zeros(img_T_R_P_h, img_T_R_P_w);
img_T_R_P_2 = zeros(img_T_R_P_h, img_T_R_P_w, nCh);
for i = 1: img_T_R_P_w
    x1_ = i - T_R_P_xy2ijShift(1);
    for j = 1: img_T_R_P_h
        y1_ = j - T_R_P_xy2ijShift(2);
        
        z1 = -Z * f * cos(theta) / (y1_ * sin(theta) - f*cos(theta));
        
        y1 = y1_ / f * z1;
        x1 = x1_ / f * z1;
        
        x0 = x1;
        y0 = y1 / cos(theta);
        img_T_R_P_toT_idx_X(j,i) = x0;
        img_T_R_P_toT_idx_Y(j,i) = y0;
    end
end

for i = 1: img_T_R_P_w
    for j = 1: img_T_R_P_h
        x0 = round(img_T_R_P_toT_idx_X(j,i) - Tij2xyShift(1));
        y0 = round(img_T_R_P_toT_idx_Y(j,i) - Tij2xyShift(2));
        x0 = min(max(x0,1), T_W);
        y0 = min(max(y0,1), T_H);
        img_T_R_P_2(j,i,:) = T(y0,x0,:);
    end
end

p = img_T_R_P_2;
return;




function template = genTemplate(varargin)
if ~exist('params')
    params = [];
end

w = 300;
h = 250;
nOfGrid = [4; 4];
gridWidth = 50;

bdH = (h-nOfGrid(2)*gridWidth)/2;
bdW = (w-nOfGrid(1)*gridWidth)/2;

template = uint8(zeros(nOfGrid(2)*gridWidth + bdH * 2, ...
    nOfGrid(1)*gridWidth + bdW * 2, 3));
for i = 1:nOfGrid(2)
    for j = mod(i,2)+1:2:nOfGrid(1)
        template(((i-1) * gridWidth+1 : i*gridWidth) + bdH, ...
            ((j-1) * gridWidth+1 : j*gridWidth) + bdW, ...
            : ) = uint8(255);
    end
end
template(bdH,3:end-2,1) = 255;
template(h-bdH+1,3:end-2,1) = 255;
template(3:end-2,bdW,1) = 255;
template(3:end-2,w-bdW+1,1) = 255;

template = double(template)/255;