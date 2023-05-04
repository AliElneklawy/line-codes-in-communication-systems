clear all;
clc;

seq = randi([0 1], 1, 10);

Rb = 1; % bit rate
Tb = 1/Rb; % bit duration
Fs = 500; % sampling frequency in Hz
f = 0:0.01:2*Rb; % frequency range
a = 1; % PSD amplitude

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

subplot(6, 1, 1);
plot(t, y);
axis([0 t(end) -2.5 2.5]);
xlabel('Time (s)');
ylabel('Amplitude');
title('Polar NRZ Modulation');
%----------------- polar NRZ ---------------%

%----------------- polar RZ ----------------%
for i = 1:length(seq)
  if seq(i) == 0
    y((i-1)*Fs+1:(i-0.5)*Fs) = -1;
    y((i-0.5)*Fs+1:i*Fs) = 0;
  else
    y((i-1)*Fs+1:(i-0.5)*Fs) = 1;
    y((i-0.5)*Fs+1:i*Fs) = 0;
  end
end

subplot(6, 1, 2);
plot(t, y);
axis([0 t(end) -2.5 2.5]);
xlabel('Time (s)');
ylabel('Amplitude');
title('Polar RZ Modulation');
%----------------- polar RZ ---------------%

%---------- Manschester coding ------------%
for i = 1:length(seq)
  if seq(i) == 1
    y((i-1)*Fs+1:(i-0.5)*Fs) = 1;
    y((i-0.5)*Fs+1:i*Fs) = -1;
  else
    y((i-1)*Fs+1:(i-0.5)*Fs) = -1;
    y((i-0.5)*Fs+1:i*Fs) = 1;
  end
end

subplot(6, 1, 3);
plot(t, y);
axis([0 t(end) -2.5 2.5]);
xlabel('Time (s)');
ylabel('Amplitude');
title('Manschester coding Modulation');
%--------- Manschester coding --------%

%----------------- AMI ---------------%
last_bit = 1; % will be used to check if the last 1-bit represented as 1 or -1
for i = 1:length(seq)
    if seq(i) == 0
        y((i-1)*Fs+1:i*Fs) = 0;
    else
        y((i-1)*Fs+1:i*Fs) = last_bit;
        last_bit = -1 * last_bit;
    end
end

subplot(6, 1, 4);
plot(t, y);
axis([0 t(end) -2.5 2.5]);
xlabel('Time (s)');
ylabel('Amplitude');
title('AMI Modulation');
%----------------- AMI ---------------%

%------- 3-level transmission --------%
last_bit = 0;
alter = [1 0 -1 0];
j = 1;
for i = 1:length(seq)
    if seq(i) == 0
        y((i-1)*Fs+1:i*Fs) = last_bit;
    else
         y((i-1)*Fs+1:i*Fs) = alter(j);
         last_bit = alter(j);
         j = j + 1;
         if j > length(alter)
             j = 1;
         end
    end
end

subplot(6, 1, 5);
plot(t, y);
axis([0 t(end) -4 2]);
xlabel('Time (s)');
ylabel('Amplitude');
title('3-Level Transmission Modulation');
%------- 3-level transmission --------%

%--------------- NRZ inverted --------------%
 last_bit = 1 ;
 for i = 1:length(seq)
    if seq(i) == 0
        y((i-1)*Fs+1:i*Fs) = last_bit;
    else
        y((i-1)*Fs+1:i*Fs) = -1 * last_bit;
        last_bit = -1 * last_bit;
    end
 end
subplot(6, 1, 6);
plot(t, y);
axis([0 t(end) -2.5 2.5]);
xlabel('Time (s)');
ylabel('Amplitude');
title('NRZ inverted Modulation');
%--------------- NRZ inverted --------------%

% Add super-title
seq_str = num2str(seq);
main_title = sprintf('Sequence: %s', seq_str);
sgtitle(main_title);

%----------------------------------
%           PSD plots
%----------------------------------
arg = f * Tb;
P_NRZ = (a^2) * Tb * sinc(arg) .* sinc(arg); 
P_RZ = (a^2 / 2) * ((sinc(arg/2)).* (sinc(arg/2)));
P_MAN = a^2 * Tb * (sinc(arg/2)).^2 .* (sin(pi*arg/2)).^2;
P_AMI = (a^2 / (2*Tb)) * (sinc(2 * pi * arg / Tb)).^2;
P_3LT = (a^2/Tb) * (0.5 * ((sin(pi*arg)./(pi*arg)).^2) + ...
    ((1./(pi*arg)).^2) .* ((sin(2*pi*arg)./2).^2));
    
figure(2)
hold on
plot(f,P_NRZ,'r')
plot(f,P_RZ,'g')
plot(f,P_MAN,'b')
plot(f,P_AMI,'k')
plot(f, P_3LT, 'm')
grid on
box on
xlabel('Frequency')
ylabel('Power Spectral Density')
title('PSD for the Line Codes')
legend('PSD for Polar NRZ Signal','PSD for polar RZ Signal',...
    'PSD for Manchester Signal', 'PSD for AMI',...
    'PSD for 3-level transmission');
