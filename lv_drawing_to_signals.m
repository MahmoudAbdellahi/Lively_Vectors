function signals = lv_drawing_to_signals(simulated_len)
% reads drawing.png and converts it to signals with length=simulated_len

drawing = imread('drawing.png');
bw = rgb2gray(drawing);

bw = bw<255;
L = bwlabel(bw);
signals=nan(max(L(:)), simulated_len);
for i=1:max(L(:))
    [r, c] = find( flip(L==i,1) );
    [~,ir] = unique(c,'rows','legacy');
    r = r(ir);
    r = imresize(r',[1 simulated_len]); r = r-mean(r);
    signals(i,:) = r;
end

end