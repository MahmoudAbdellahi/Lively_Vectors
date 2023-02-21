function signal = lv_imgtosig
% converts image of a signal to actual signal byt thersholding the pixel
% values ... the length of the signal is 1000samples and the limits are
% always 0 to 1
[file,path] = uigetfile;

A = double(rgb2gray(imread([path file])));
A = flip(A);

[y,x] = find(A>200);

new_x = interp1([min(x) max(x)],[0 1],x);
new_y = interp1([min(y) max(y)],[0 1],y);

% shape trimmer to trim duplications and make it more signal like (because thick line will give many pixels)
y2=[];
for i=1:length(new_x)
    y2 = [y2 ; median(new_y( new_x==new_x(i) ))];
end
x2= (new_x);

xy = unique([x2 y2],'rows');


x_query = linspace(0,1,1000);
new_y = interp1(xy(:,1),xy(:,2),x_query);

signal = new_y;    
  
plot(x_query,signal);

end
