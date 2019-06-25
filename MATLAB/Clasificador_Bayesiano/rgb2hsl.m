function HSL = rgb2hsl(RGB)
[lv, lu] = size(RGB(:,:,1));
RGB = im2double(RGB);
HSL = zeros(lv,lu,3);
for i1 = 1:lv
    for i2 =1:lu
        [H,S,L] = HSLPIX(RGB(i1,i2,:));
        HSL(i1,i2,1) = H;
        HSL(i1,i2,2) = S;
        HSL(i1,i2,3) = L;
    end
end 
end