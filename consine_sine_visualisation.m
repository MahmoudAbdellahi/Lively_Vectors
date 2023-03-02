a = 1; % amplitude
fs=100; % sampling rate
t = 0:1/fs:2; % time points
f = 2; % frequency (cycles per second)
phi = 0; % phase shift

signal_sine = a*sin(2*pi*f*t + phi);
signal_cosine = a*cos(2*pi*f*t + phi);

plot(t, signal_cosine) % time, amplitude
%% static points
figure,
for i=1:round(length(t)/8):length(t)
    subplot(2,2,1) 
    [u,v] = pol2cart(2*pi*f*t(i), 1);
    c = compass(u,v);  
    hold on, plot(signal_cosine(i),0,'ob'), hold off, % cosine projection
    hold on, plot(0,signal_sine(i),'ob'), hold off, % sine projection
    c1 = c(1); c1.LineWidth = 2; c1.Color = 'r';  axis square;
    set(findall(gcf, 'String', '0'),'String', 'cosine'); 
    set(findall(gcf, 'String', '90'),'String', 'sine');  drawnow;

    subplot(2,2,2)
    hold on, plot(t(i), signal_sine(i),'o','color','k'), xlim([t(1) t(end)]), ylim([-a a]), xlabel('sine wave projection'), ylabel('amplitude'), axis square;

    subplot(2,2,3)
    hold on, plot(t(i), signal_cosine(i),'o','color','k'), xlim([t(1) t(end)]), ylim([-a a]), xlabel('cosine wave projection'), ylabel('amplitude'), axis square;
end


%% all points
flag=0;
for i=1:length(t)
    subplot(2,2,1) 
    [u,v] = pol2cart(2*pi*f*t(i), 1);
    c = compass(u,v);  
    hold on, plot(signal_cosine(i),0,'ob'), hold off, % cosine projection
    hold on, plot(0,signal_sine(i),'ob'), hold off, % sine projection
    c1 = c(1); c1.LineWidth = 2; c1.Color = 'r';  axis square;
    set(findall(gcf, 'String', '0'),'String', 'cosine'); 
    set(findall(gcf, 'String', '90'),'String', 'sine');  drawnow;
    if 2*pi*f*t(i)>=1.0472 && flag==0, cosi=signal_cosine(i), si=signal_sine(i), flag=1; end % 60 degrees

    subplot(2,2,2)
    hold on, plot(t(i), signal_sine(i),'o','color','k'), xlim([t(1) t(end)]), ylim([-a a]), xlabel('sine wave projection'), ylabel('amplitude'), axis square;

    subplot(2,2,3)
    hold on, plot(t(i), signal_cosine(i),'o','color','k'), xlim([t(1) t(end)]), ylim([-a a]), xlabel('cosine wave projection'), ylabel('amplitude'), axis square;
    pause(0.1)
end

%% 3d with angles as third dimension
figure
ax = nexttile; rot=0;
rotx = 10; roty=10; rotz=20;
view(ax,[rotx roty rotz]);
for i=1:length(t)
    hold on, 
    plot3(ax, signal_cosine(i), signal_sine(i), t(i),'o','color','k'), xlim([-a a]), ylim([-a a]), zlim([t(1) t(end)]);
    xlabel('cosine'), ylabel('sine'), zlabel('temporal progression'), axis square;
    if rotx<=100 % we want 90 and we started at 10
        rotx = rotx+2; roty = roty-0.5; if roty<=0; roty=0; end
        rotz = rotz-0.5; if rotz<=0; rotz=0; end
        view(ax,[rotx roty rotz]); 
    end

    pause(0.1)
end

 




