function out = plotImgPoint(I, pos, wh)

[h, w, ch] = size(I);
tmp = I;
I = zeros(h+2,w+2,ch);
I(2:end-1,2:end-1,:) = tmp;
h = h+2;
w = w+2;
pos = pos+1;

ow = max(pos(1,:));
oh = max(pos(2,:));
out = zeros(wh(2),wh(1),ch);
out_tmp = zeros(wh(2),wh(1));

x_b = min(max(floor(pos(1,:)),1),w);
x_u = min(max(ceil(pos(1,:)),1),w);
x_w_b = 1 - min(abs(pos(1,:) - x_b), 1);
x_w_u = 1 - x_w_b;

y_b = min(max(floor(pos(2,:)),1),h);
y_u = min(max(ceil(pos(2,:)),1),h);
y_w_b = 1 - min(abs(pos(2,:) - y_b), 1);
y_w_u = 1 - y_w_b;

for i = 1:ch
    in_tmp = I(:,:,i);
    
    tmp = floor(((1:length(x_b))-1)/length(x_b)) * w*h;
    idx_ybxb = (x_b-1)*h+y_b + tmp;
    idx_ybxu = (x_u-1)*h+y_b + tmp;
    idx_yuxb = (x_b-1)*h+y_u + tmp;
    idx_yuxu = (x_u-1)*h+y_u + tmp;
    w_ybxb = x_w_b.*y_w_b;
    w_ybxu = x_w_u.*y_w_b;
    w_yuxb = x_w_b.*y_w_u;
    w_yuxu = x_w_u.*y_w_u;
    
    out_tmp(:) = in_tmp(idx_ybxb) .* w_ybxb + ...
        in_tmp(idx_ybxu) .* w_ybxu + ...
        in_tmp(idx_yuxb) .* w_yuxb + ...
        in_tmp(idx_yuxu) .* w_yuxu;
    
    out(:,:,i) = out_tmp;
end





