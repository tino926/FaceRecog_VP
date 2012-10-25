function [p, template] = genPerspective(params)

template = genTemplate();
[sTemplateH, sTemplateW,nChannel] = size(template);


Z = 100;
f = 100;
theta = -pi/50;

zeroCenterToTemplate = [(sTemplateW+1)/2; (sTemplateH+1)/2];
xOnTemplate = repmat(1:sTemplateW, sTemplateH, 1);
xOnTemplate = xOnTemplate(:)';
yOnTemplate = repmat(1:sTemplateH, 1, sTemplateW);

x3d = xOnTemplate - zeroCenterToTemplate(1);
y3d = (yOnTemplate - zeroCenterToTemplate(2)) * cos(theta)+ sin(theta) * Z;
z3d = cos(theta) * Z - (yOnTemplate - zeroCenterToTemplate(2))* sin(theta);

xProjected3d = x3d./z3d*f;
yProjected3d = y3d./z3d*f;

zeroCenterToProjected3d = [1;1] - [min(xProjected3d); min(yProjected3d)];


xOnProjected = round(xProjected3d+zeroCenterToProjected3d(1));
yOnProjected = round(yProjected3d+zeroCenterToProjected3d(2));
sProjectedW = max(xOnProjected);
sProjectedH = max(yOnProjected);

% for i = 1:nChannel
%     tmp1 = template(:,:,i);
%     tmp2 = p(:,:,i);
%     tmp2( (xOnProjected-1)*sProjectedH + yOnProjected) = tmp1(:);
%     p(:,:,i) = tmp2;
% end

% reverse
p = uint8(zeros(sProjectedH, sProjectedW, nChannel));
xOnProjected_ = repmat(1:sProjectedW, sProjectedH, 1);
xOnProjected_ = xOnProjected_(:)';
yOnProjected_ = repmat(1:sProjectedH, 1, sProjectedW);

y3d


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
template(bdH,:,1) = 255;
template(h-bdH+1,:,1) = 255;
template(:,bdW,1) = 255;
template(:,w-bdW+1,1) = 255;