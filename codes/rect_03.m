% clear all;
% x:horizontal axis, y:vertical axis,

tic
clear all
DataPath = 'D:\WorkingData\FaceRecog_VP';

l = [253, 132, 259, 210; 373 18 371 126; 541 50 522 197];
img = double(imread(fullfile(DataPath,'.\2_640x480\C1_000000328.bmp')))/255;
% l = [320, 373, 338, 481; 34 560 80 665; 1149 279 1079 498];
% img = double(imread('pic.JPG'))/255;
% l = [111, 95, 172, 412; 383 150 379 274; 482 189 467 306];
% img = double(imread('Picture 15.JPG'))/255;
% l = [227, 133, 251, 358; 513 134 497 346];
% img = double(imread('sync01.png'))/255;
% l = [80, 197, 122, 368; 207 212 220 312; 598 137 550 330];
% img = double(imread('sync02.png'))/255;
% [img, l] = genPerspective01();
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

% parameter of central separation line
tmp = (IC+VP)/2;

dPVIC = VP-IC;
dPVIC_normal = normc(dPVIC);

% mRotateInv * VP and mRotateInv * IC should on the same vertical line
mRotateInv = [dPVIC_normal(2) -dPVIC_normal(1);...
    dPVIC_normal(1) dPVIC_normal(2)];
mRotate = [dPVIC_normal(2) dPVIC_normal(1);...
    -dPVIC_normal(1) dPVIC_normal(2)];

tmp = [1 w 1 w; 1 1 h h] - repmat(IC, [1 4]);
tmp = mRotateInv*tmp;

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

img_R_I = zeros(img_R_I_h, img_R_I_w, 3);



% <-- rectify rotation
x = repmat( ((1:img_R_I_w) - img_R_I_xy2ijShift(1)), ...
    img_R_I_h, 1) + IC(1);
x = x(:)';
y = repmat( ((1:img_R_I_h) - img_R_I_xy2ijShift(2)), ...
    1, img_R_I_w) + IC(2);

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

img_R_I(:) = img(idx_ybxb) .* w_ybxb + ...
    img(idx_ybxu) .* w_ybxu + ...
    img(idx_yuxb) .* w_yuxb + ...
    img(idx_yuxu) .* w_yuxu;
imshow(img_R_I);
% --> <End> rectify rotation


% img_rotateInvLine = img_R_I;
% 
% for i = 10:30:size(img_rotateInvLine,2)-1
%     a_tmp = (i - img_R_I_VP(1)) / (1 - img_R_I_VP(2));
%     b_tmp = i - a_tmp;
%     
%     for j = 1:size(img_rotateInvLine,1);
%         img_rotateInvLine( j, round(a_tmp*j + b_tmp), 1) = 255;
%         img_rotateInvLine( j, round(a_tmp*j + b_tmp), 2) = 0;
%         img_rotateInvLine( j, round(a_tmp*j + b_tmp), 3) = 0;
%     end
% end




tmpx = zeros(img_R_I_h, img_R_I_w);
tmpy = zeros(img_R_I_h, img_R_I_w);

% now suppose the center is on img_rotateInv's top-middle
zeroCenter = [0; 0];
shiftToZeroCenter = [zeroCenter(1)-img_R_I_IC(1); 0];
IC_zeroCenter = img_R_I_h + shiftToZeroCenter;
VP_zeroCenter = img_R_I_VP + shiftToZeroCenter;

img_rotateInvLine = img_R_I;

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
figure, imshow(img_rotateInvLine);




% assum Z and f, from f and VP, we can get theta
Z = 600;
f = 600;

img_R_I_IC_xy = img_R_I_IC - img_R_I_IC;
img_R_I_VP_xy = img_R_I_VP - img_R_I_IC;
VP_dY = img_R_I_VP_xy(2);
theta = atan(f/VP_dY);
theta_ground = theta - pi/2;


% map from img_R_I to img_rect
% for i = 1:img_R_I_w
%     for j = 1:img_R_I_h
%         
%     end
% end


















img_rect_w = 501;
img_rect_h = 401;
img_rect_ij2xyShift = -[(img_rect_w+1)/2; (img_rect_h+1)/2];

img_rect = zeros(img_rect_h, img_rect_w,3);
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

img_rect_R_P_idx_Xij = round(img_rect_R_P_idx_X+img_R_I_IC(1));
img_rect_R_P_idx_Yij = round(img_rect_R_P_idx_Y+img_R_I_IC(2));

img_rect_R_P_idx_Xij = min( img_R_I_w, max(1,img_rect_R_P_idx_Xij));
img_rect_R_P_idx_Yij = min( img_R_I_h, max(1,img_rect_R_P_idx_Yij));

for i = 1: img_rect_w
    for j = 1: img_rect_h
        x1_ = round(img_rect_R_P_idx_Xij(j,i));
        y1_ = round(img_rect_R_P_idx_Yij(j,i));
        img_rect(j,i,:) = img_rotateInvLine(y1_,x1_,:);
    end
end

figure, imshow(img_rect)

















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


















