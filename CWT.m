Fs = 1/(time(2)-time(1));
wp=[3,120];
sample_length = 746; %number of data points per CA pulse = pulse width / sampling interval: 746 744

% [cfs,f] = cwt(dvdt,Fs,'WaveletPrameters', wp,'VoicesPerOctave', 6, 'FrequencyLimits', [1.5 Fs/2]);
% [ccs,f] = cwt(didt,Fs','WaveletParameters', wp, 'VoicesPerOctave', 6, 'FrequencyLimits', [1.5 Fs/2]);

[cfs,f] = cwt(potential./10,Fs,'WaveletPrameters', wp,'VoicesPerOctave', 6, 'FrequencyLimits', [5 Fs]);
[ccs,f] = cwt(current,Fs','WaveletParameters', wp, 'VoicesPerOctave', 6, 'FrequencyLimits', [5 Fs]);

% [cfs,f] = cwt(potential./10,'cmor1-.45',Fs);
% [ccs,f] = cwt(current,'cmor1-.45',Fs');

% [cfs,f] = cwt(potential./10,Fs, 'amor','VoicesPerOctave', 4, 'FrequencyLimits', [1.2 Fs/2]);
% [ccs,f] = cwt(current,Fs','amor','VoicesPerOctave', 4, 'FrequencyLimits', [1.2 Fs/2]);

% [cfs,f] = cwt(potential./10,Fs, 'bump');
% [ccs,f] = cwt(current,Fs','bump');  'cmor1-.45'


imp_cwt=cfs./ccs;

%% scalogram
 figure(1);
 clf;
 helperHyperbolicChirpPlot(cfs,f,time(1:end)); %-1
 figure(2);
 clf;
 helperHyperbolicChirpPlot(ccs,f,time(1:end));

%% CA plot
figure(4);
clf;
plot(time, current);

%% EIS comparison
figure(5);
ca_no = 3;
ca_no_2 = 71;
plot(real(imp_cwt(1:end,sample_length*2*ca_no)), -imag(imp_cwt(1:end,sample_length*2*ca_no)),'ko-'); 
hold on;
plot(real(imp_cwt(1:end,sample_length*2*ca_no_2)), -imag(imp_cwt(1:end,sample_length*2*ca_no_2)),'ro-');
% plot(real(imp_cwt(1:end,sample_length*2*42)), -imag(imp_cwt(1:end,sample_length*2*42)),'bo');
% plot(real(imp_cwt(1:end,sample_length*2*66)), -imag(imp_cwt(1:end,sample_length*2*66)),'mo');
hold off;

%% scalogram of phase angle and current overlap
angle_imp = angle(imp_cwt)*180/pi;
imp_cwt_flip = flip(imp_cwt,1);
real_imp = real(imp_cwt)./1e6;
idxtodel = (current < 0) | (current > 4e-8);
cu_del= current;
ti_del = time;
cu_del(idxtodel) = [];
ti_del(idxtodel) = [];

figure(63)
clf;
imagesc(time, log2(f), -angle_imp);
colorbar;
set(gca,'YDir','normal');
colormap("jet");
caxis([-50 50]);
ylim([log2(5.53) log2(7.78)]);
yticks([log2(5.83) log2(6.54) log2(7.35) log2(8)]);
%styling
cb = colorbar;
%cb.Label.String = 'Z^\prime (M\Omega)';
cb.Label.String = '-Phase angle (Degree)';
cb.FontSize = 20;
cb.Label.FontSize = 24;
cb.Label.FontWeight = 'bold';
xlim([10 16]);
ax = gca;
ax.YTickLabels = 2.^yticks;
ax.FontSize = 20;
xlabel('Time (s)' , 'FontSize', 24, 'FontWeight', 'bold');
ax.YLabel.String = 'Frequency (Hz)';
ax.YLabel.FontSize = 24;
ax.YLabel.FontWeight = 'bold';
%ylabel('Frequency (Hz)' ,  'FontSize', 24, 'FontWeight', 'bold');


hold on
ratio = 0.2;
pos = 2.45; %2.7
plot(ti_del, (cu_del.*ratio*1e8)+pos ,'k','LineWidth',3); %'color',"#ffa500"
hold off

%% bode plot
 figure(7);clf;
 for i = 1:74
     subplot(2,1,1);
     semilogx(f, (log10(abs(imp_cwt(1:end,sample_length*2*i)))));
     xlabel('log frequency (Hz)','FontSize', 14, 'FontWeight', 'bold');
     ylabel('log|Z| (\Omega)','FontSize', 14, 'FontWeight', 'bold');
     hold on;
     subplot(2,1,2);
     semilogx(f, -angle(imp_cwt(1:end,sample_length*i))*180/pi);
     xlabel('log frequency (Hz)','FontSize', 14, 'FontWeight', 'bold');
     ylabel('phase (deg)','FontSize', 14, 'FontWeight', 'bold');
     hold on;
 end
 hold off;
 
%% 3d nyquist plot
real_imp_cwt = real(imp_cwt);
imag_imp_cwt = -imag(imp_cwt);

idxtodel1 = real_imp_cwt < 0 | real_imp_cwt > 7e6;
idxtodel2 = imag_imp_cwt > 3e6 | imag_imp_cwt < 0 ;

real_del_cwt = real_imp_cwt./1e6;
imag_del_cwt = imag_imp_cwt./1e6;

real_del_cwt(idxtodel1) = NaN;   %[];
imag_del_cwt(idxtodel2) = NaN;   %[];
figure(8);
arr = zeros(length(imp_cwt(1:end,1)));

for i = 1:(74)
    %plot3(arr+i, real(imp_cwt(1:end,sample_length*2*i)), -imag(imp_cwt(1:end,sample_length*2*i)),'k.'); % N = pulse width / sampling interval
    plot3(arr+i, real_del_cwt(1:end,sample_length*2*i), imag_del_cwt(1:end,sample_length*2*i),'k.');
    hold on;
end
grid on;
hold off;
axis on;
%axis([0 80 0 inf 0 inf]);
xlabel('Time (s)','FontSize', 24, 'FontWeight', 'bold');
ylabel('Z^\prime (M\Omega)','FontSize', 24, 'FontWeight', 'bold');
zlabel('-Z^\prime^\prime (M\Omega)','FontSize', 24, 'FontWeight', 'bold');
%xticks([0 20 40 60 80]);
yticks([0 1.5 3 4.5  6]);
zticks([0 1 2 3 4]);
ax = gca;
ax.FontSize = 16;
axis([0 inf 0 6 0 4]);

%% 2d nyquist plot
figure(9);
for i = 1:74
    plot(real_del_cwt(1:end,sample_length*2*i), imag_del_cwt(1:end,sample_length*2*i),'k.'); % N = pulse width / sampling interval
    hold on;
end
hold off;
axis equal;


%% save data
num_nyquist = 4;
reZ_cwt = real(imp_cwt(1:end,sample_length*2*num_nyquist));
imZ_cwt = imag(imp_cwt(1:end,sample_length*2*num_nyquist));
f_cwt = flip(f);
data_Z = [f reZ_cwt imZ_cwt];

%% taking imp data at specific frequency
for i = 1:74 %74
    time_plot(i) = time(sample_length*2*i)';
    abs_plot(i) = abs(imp_cwt(end,sample_length*2*i));
    real_plot(i) = real(imp_cwt(end,sample_length*2*i)); %real_plot(i) = real(imp_cwt(end,sample_length*2*i)); real_plot(i) = real(imp_cwt(end-2,sample_length*i));
    %mag_plot(i) = imag(imp_cwt(1,sample_length*2*i));
end
time_plot_t = time_plot';
abs_plot_t =  abs_plot';
real_plot_t = real_plot';
figure(14)
plot(time_plot, real_plot, '-ko');

%% whole data set
resol = time(2)-time(1);
time_re_del = (1:1:length(real_del_cwt))'.*resol;
