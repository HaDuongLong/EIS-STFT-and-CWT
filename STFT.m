fs=1/(time(2)-time(1));

opts = {'Window',kaiser(1002,10),'FrequencyRange',"centered"}; %,'FrequencyRange'

opts1= {'Window',hamming(1492,'periodic'),'OverlapLength', 0 ,'FrequencyRange',"onesided"};

sample_length = 746;

[Sc,Fc,Tc] = stft(current,fs,opts1{:});
[Sp,Fp,Tp] = stft(potential./10,fs,opts1{:});
[Si,Fi,Ti] = stft(potential./10,fs,opts1{:});


for i=1:length(Sc(1,:))
    Imp_tt(:,i)=Sp(:,i) ./ Sc(:,i);
end
Imp_tt_t = Imp_tt';
Fp_pos = Fp(Fp>0);
Fp_pos_f = flip(Fp_pos);

%% Phase Spectrogram and current overlay
figure(16)
fp_log2 = log2(Fp(5:end));
angle_imp_tt = angle(Imp_tt(5:end,1:end))*180/pi;
imagesc(Tp, Fp, -angle_imp_tt); %real(Imp_tt)./1e6
colormap("jet");
colorbar;
cb = colorbar;
cb.Label.String = '-Phase angle (Degree)'; %'Z^\prime (M\Omega)'
cb.FontSize = 20;
cb.Label.FontSize = 24;
cb.Label.FontWeight = 'bold';
ax = gca;
ax.FontSize = 20;
axis tight
xlabel('Time (s)','FontSize', 24, 'FontWeight', 'bold');
ylabel('Frequency (Hz)','FontSize', 24, 'FontWeight', 'bold');
xlim([10 16]);
set(gca,'YDir','normal')
ylim([4.5 8.5]);
yticks([ 5 6 7 8]);
caxis([-50 50]); %phase
hold on
idxtodel = (current < 0) | (current > 4e-8);
cu_del= current;
ti_del = time;
cu_del(idxtodel) = [];
ti_del(idxtodel) = [];
ratio = 1.6;
pos = 4.5; %2.7
plot(ti_del, (cu_del.*ratio*1e8)+pos ,'k','LineWidth',3); %'color',"#ffa500"
% axis([0 75 0 inf]);
hold off

%% 3D nyquist plot
figure(3);
%delete out of range value
real_imp_tt = real(Imp_tt);
imag_imp_tt = -imag(Imp_tt);
idxtodel1 = real_imp_tt < 0 | real_imp_tt > 7e6;
idxtodel2 = imag_imp_tt > 4e6 | imag_imp_tt < 0 ;
real_del = real_imp_tt;
imag_del = imag_imp_tt;
real_del(idxtodel1) = NaN;   %[];
imag_del(idxtodel2) = NaN;   %[];

plot3(Tp, real_del(4:5:400,:)./1e6, imag_del(4:5:400,:)./1e6 ,'k.');
grid on;
axis([0 80 0 6 0 4]);
xlabel('Time (s)','FontSize', 24, 'FontWeight', 'bold');
ylabel('Z^\prime (M\Omega)','FontSize', 24, 'FontWeight', 'bold');
zlabel('-Z^\prime^\prime (M\Omega)','FontSize', 24, 'FontWeight', 'bold');
xticks([0 20 40 60 80]);
yticks([0 1.5 3 4.5  6]);
zticks([0 1 2 3 4]);
ax = gca;
ax.FontSize = 16;

%% plot imps at specific freq index 6 = 5Hz
for i = 1:75
    abs_plot(i) = abs(Imp_tt(6,i));
    real_plot(i) = real(Imp_tt(24,i));
end
figure(9)
plot(Ti(2:75), real_plot(2:end),'-ko');

real_plot_t =  real_plot';

%% EIS comparison
figure(6);
ca_no = 11;
ca_no_2 = 14;
plot(real_del(2:2:350,ca_no)./1e6, imag_del(2:2:350,ca_no)./1e6 ,'ko-');
hold on
plot(real_del(2:2:350,ca_no_2)./1e6, imag_del(2:2:350,ca_no_2)./1e6 ,'ro-');
hold off
