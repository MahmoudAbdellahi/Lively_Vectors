%% Fourier transform of a real valued wave
clear all
% signal = lv_drawing_to_signals(200); % real valued signal won't leave 
% the real xaxis in: real, imaginary coordinates
fs = 100; 
t = 0:1/fs:2-1/fs; 
signal = 4 * cos(2*pi*2*t); 
figure,
plot(t, signal); % 100 samples per second



N = length(signal); 
complex_sine_wave_t = (0:N-1)/N; % normalised time 
fr = 0:N-1; % to always have frequencies that are independent of the unit..

nyquist = fs/2; 
pts_hz = linspace(0,nyquist,(N/2)+1); % cannot go higher to prevent aliasing

for i=1:length(fr)
    complex_sine_wave = exp(-2*pi*1i*fr(i)*complex_sine_wave_t); % stop it to see the magnitude
    coeff(i) = complex_sine_wave*signal';
end
coeff = coeff./N; 
m = abs(coeff);
m(2:length(pts_hz)) = (m(2:length(pts_hz)))*2; % doubling the positive frequencies
plot(pts_hz, m(1:length(pts_hz))); 
% figure, % using matlab's fft
% temp = fft(signal);
% plot(pts_hz, (abs(temp(1:length(pts_hz)))*2)/N ); 

% the inverse fourier transform
reconstructed = zeros(1,length(signal));
for i=1:length(fr)
    reconstructed = reconstructed + (coeff(i) * exp(2*pi*1i*fr(i)*complex_sine_wave_t));
end
plot(t, real(reconstructed)) 
%% plotting the reconstructed signals with the coefficients/vectors that build it
figure, 
for j=1:length(complex_sine_wave_t) % all frequencies and specific time
    temp = coeff .* exp(2*pi*1i*fr*complex_sine_wave_t(j));  

    subplot(121), % we want to see the rotating vectors that's why we will use the yaxis for amplitude (sum of vectors) 
    % and the xaxis to extend the vectors and see them, we use as many
    % points as we have for the signal to draw the vectors..
    v = [repmat(j,length(temp),1) real(temp)']; origin=[j 0];
    plot_added_vectors(origin, v), title('rotating vectors')
    xlim([length(t) length(t)*2]), ylim([-300 300]); % xaxis is from length(t) for visualisation but we can start from 0 to see all vectors

    subplot(122),  drawnow,
    plot(j, real(reconstructed(j)),'bo'); hold on, title('constructed signal')
    xlim([0 length(t)]), ylim([-300 300]);
end



%% let's try a 3d version with complex signal 
% now we visualise real, imaginary, time progression

% cosine + i sine and then we see the vectors in 3d
% because we will plot real imag and then time pts .. 
fs=100; 
t=0:1/fs:1; freq = 2; a=3;
signal = a*cos(2*pi*freq.*t) + 1i*a*sin(2*pi*freq.*t); % complex wave
nyquist = fs/2; N = length(signal); pts_hz = linspace(0,nyquist,(N/2)+1);
plot(t,real(signal), t,imag(signal))
complex_sine_wave_t = (0:N-1)/N;

% fft
figure, % using matlab's fft
coeff = fft(signal);
plot(pts_hz, (abs(coeff(1:length(pts_hz)))*2)/N ); 
coeff = coeff./N;

% ifft
figure, % using matlab's fft
reconstructed = ifft(coeff.*N);
plot(t,real(reconstructed), t,imag(reconstructed))

% plotting the reconstructed signals with the coeffients/vectors that build it
figure, 
for j=1:length(complex_sine_wave_t) % all frequencies and specific time
    temp = coeff .* exp(2*pi*1i*freq*complex_sine_wave_t(j));  
    
    v = [real(temp)' imag(temp)'];
    v = cumsum(v,1); h=[]; origin=[0 0 0];
    for i=1:size(v,1)
        hold on,
        h = [h plot3([origin(1) ,v(i,1)], [origin(2) v(i,2)], [origin(3) ,j], '-')];
        origin = [v(i,:) j];
    end
    fprintf(['finished plotting: ' num2str(i) ' vectors. \n'])
    pause(0); set(h,'Visible','off'); % call time

    tip = sum(temp);
    hold on,
    plot3(real(tip), imag(tip), j, 'o'); view(-123,47),
    xlim([-3 3]), ylim([-3 3]); zlim([1 length(complex_sine_wave_t)]);
    xlabel('real'); ylabel('imaginary'); zlabel('time progression'); drawnow
end



 




%% helping function
function plot_added_vectors(origin, v)
% takes v (vectors by coord.) and draws all vectors visually added together
v = cumsum(v,1); h=[];
v(:,1)= v(1,1) + (0:size(v,1)-1)'; % just for scaling the xaxis so we draw the vectors in a number of points that matches the length of the signal
for i=1:size(v,1)
    hold on,
    h = [h plot([origin(1) ,v(i,1)], [origin(2) v(i,2)], 'b')]; 
    origin = v(i,:);
end
pause(0.0); set(h,'Visible','off'); % call time
hold on,
plot([origin(1) ,v(i,1)], [origin(2) v(i,2)], 'ro')

end



