clear all;
clc;

seq = randi([0 1], 1, 10);

Tb = 1; % bit duration
Fs = 100; % sampling frequency in Hz

t = linspace(0, Tb*length(seq), length(seq)*Fs+1);
t = t(1:end-1); % remove the last point to make sure that t and y have
                % the same length

%----------------- polar NRZ ---------------%
for i = 1:length(seq)
    if seq(i) == 0
        y((i-1)*Fs+1:i*Fs) = -1;
    else
        y((i-1)*Fs+1:i*Fs) = 1;
    end
end

% Plot the modulated signal
subplot(6, 1, 1);
plot(t, y);
axis([0 t(end) -2.5 2.5]);
xlabel('Time (s)');
ylabel('Amplitude');
title('Polar NRZ Modulation');
%----------------- polar NRZ ---------------%

%----------------- polar RZ ---------------%
for i = 1:length(seq)
  if seq(i) == 1
    y((i-1)*Fs+1:(i-0.5)*Fs) = 1;
    y((i-0.5)*Fs+1:i*Fs) = 0;
  else
    y((i-1)*Fs+1:(i-0.5)*Fs) = -1;
    y((i-0.5)*Fs+1:i*Fs) = 0;
  end
end


% Plot the modulated signal
subplot(6, 1, 2);
plot(t, y);
axis([0 t(end) -2.5 2.5]);
xlabel('Time (s)');
ylabel('Amplitude');
title('Polar RZ Modulation');
%----------------- polar RZ ---------------%


