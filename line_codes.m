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


% Add super-title
seq_str = num2str(seq);
main_title = sprintf('Sequence: %s', seq_str);
sgtitle(main_title);
