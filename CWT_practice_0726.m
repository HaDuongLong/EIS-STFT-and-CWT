Fs = 1/(time(2)-time(1));
wp=[3,120];
dt = diff(time);
di = diff(current);
dv = diff(potential./10);
sample_length = 746; %number of data points per CA pulse = pulse width / sampling interval: 746 744
dvdt = dv./dt;
didt = di./dt;

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
% figure(11)
% cwt(current,Fs);
% colormap("jet");
%caxis([0 2]);
% figure(12)
% cwt(potential,Fs);
% colormap("jet");
%%
% figure(1);
% clf;
% helperHyperbolicChirpPlot(cfs,f,time(1:end)); %-1
% figure(2);
% clf;
% helperHyperbolicChirpPlot(ccs,f,time(1:end));
%%
% figure(3);
% clf;
% imagesc(time,f,abs(ccs));
%caxis([3e6 10e7]);
%axis([0 25 100 1500 1e6 10e7]);   
%helperHyperbolicChirpPlot(log(abs(imp)),f,time);
%%
figure(4);
clf;
plot(time, current);
%%
% figure(10)
% clf;
% plot(potential.*10, current);
%%
figure(5);
ca_no = 3;
ca_no_2 = 71;
%plot(imp(:,10000),'k.');
plot(real(imp_cwt(1:end,sample_length*2*ca_no)), -imag(imp_cwt(1:end,sample_length*2*ca_no)),'ko-'); % N = pulse width / sampling interval
hold on;
plot(real(imp_cwt(1:end,sample_length*2*ca_no_2)), -imag(imp_cwt(1:end,sample_length*2*ca_no_2)),'ro-');
% plot(real(imp_cwt(1:end,sample_length*2*42)), -imag(imp_cwt(1:end,sample_length*2*42)),'bo');
% plot(real(imp_cwt(1:end,sample_length*2*66)), -imag(imp_cwt(1:end,sample_length*2*66)),'mo');
%plot(real(imp_cwt(2:70,2906)), -imag(imp_cwt(2:70,2906)),'r.');
hold off;
%plot(imp(:, 0.1*Fs));
%%
angle_imp = angle(imp_cwt)*180/pi;
imp_cwt_flip = flip(imp_cwt,1);
figure(22);
clf;
%subplot(2,1,1)
% imagesc(time(1:end),f,-angle_imp);
% set(gca,'YDir','normal');
% colorbar;
% flip_angle = flipud(angle_imp);
real_imp = real(imp_cwt)./1e6;
helperHyperbolicChirpPlot1(angle_imp,f,time(1:end)); %./1e6 real_imp
colormap("jet");
%colormap(flipud(angle_imp));
%colorbar off;
ylim([log2(5.830224045931709) log2(13.394337541651836)]); % 2.56 3.7
%yticks([log2(5.830224045931709) log2(7.693026750185245) log2(10.151009654656920) log2(13.394337541651836)]);
%caxis([1.5 6.9]);
caxis([-120 0]);
%ylim([3.15 3.2]);
% xlim([10 16]);
% ax = gca;
% ax.TickDir = 'out';
% ax.TickLength = [0.01 0.01];

hold on
idxtodel = (current < 0) | (current > 4e-8);
cu_del= current;
ti_del = time;
cu_del(idxtodel) = [];
ti_del(idxtodel) = [];
%subplot(2,1,2)
ratio = 0.27;
pos = 2.7; %2.7
plot(ti_del, (cu_del.*ratio*1e8)+pos ,'b','LineWidth',3); %'color',"#ffa500"
hold off
%%
figure(63)
clf;
imagesc(time, log2(f), -angle_imp);
colorbar;
set(gca,'YDir','normal');
colormap("jet");
%axis limit
%caxis([1.5 6.9]);
caxis([-50 50]);
% ylim([log2(5.53) log2(7.78)]);
% yticks([log2(5.83) log2(6.54) log2(7.35) log2(8)]);
%styling
cb = colorbar;
%cb.Label.String = 'Z^\prime (M\Omega)';
cb.Label.String = '-Phase angle (Degree)';
cb.FontSize = 20;
cb.Label.FontSize = 24;
cb.Label.FontWeight = 'bold';
%xlim([10 16]);
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
%% 
% figure(62)
% clf;
% imagesc(time(1:end),f,(real(imp_cwt)));
% colormap("jet");
% caxis([1e5 0.4e7]);
%% bode plot
% figure(6);clf;
% ca_no1 = 2;
% subplot(2,1,1);
% semilogx(f, (log10(abs(imp_cwt(1:end,sample_length*2*ca_no1)))),'k');
% xlabel('log frequency (Hz)','FontSize', 14, 'FontWeight', 'bold');
% ylabel('log|Z| (\Omega)','FontSize', 14, 'FontWeight', 'bold');
% subplot(2,1,2);
% semilogx(f, -angle(imp_cwt(1:end,sample_length*2*ca_no1))*180/pi,'k');
% xlabel('log frequency (Hz)','FontSize', 14, 'FontWeight', 'bold');
% ylabel('phase (deg)','FontSize', 14, 'FontWeight', 'bold');

%% bode plot
% figure(7);clf;
% for i = 1:20
%     subplot(2,1,1);
%     semilogx(f, (log10(abs(imp_cwt(1:end,sample_length*2*i)))));
%     xlabel('log frequency (Hz)','FontSize', 14, 'FontWeight', 'bold');
%     ylabel('log|Z| (\Omega)','FontSize', 14, 'FontWeight', 'bold');
%     hold on;
%     subplot(2,1,2);
%     semilogx(f, -angle(imp_cwt(1:end,sample_length*i))*180/pi);
%     xlabel('log frequency (Hz)','FontSize', 14, 'FontWeight', 'bold');
%     ylabel('phase (deg)','FontSize', 14, 'FontWeight', 'bold');
%     hold on;
% end
% hold off;
 
%% 3d
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
% for i = 1:2:(74*2)
%     %plot3(arr+i, real(imp_cwt(1:end,sample_length*2*i)), -imag(imp_cwt(1:end,sample_length*2*i)),'k.'); % N = pulse width / sampling interval
%     plot3(arr+i, real_del_cwt(1:end,sample_length*i), imag_del_cwt(1:end,sample_length*i),'k.');
%     hold on;
% end
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
%annotation(figure(8),'textbox','String',{'A'},'FontSize', 24, 'FontWeight', 'bold','LineStyle','none');
axis([0 inf 0 6 0 4]);

%% 2d
figure(9);
for i = 1:74
    plot(real_del_cwt(1:end,sample_length*2*i), imag_del_cwt(1:end,sample_length*2*i),'k.'); % N = pulse width / sampling interval
    hold on;
end
hold off;
axis equal;
%axis([0 inf 0 inf]);

%% save data
num_nyquist = 4;
reZ_cwt = real(imp_cwt(1:end,sample_length*2*num_nyquist));
imZ_cwt = imag(imp_cwt(1:end,sample_length*2*num_nyquist));
f_cwt = flip(f);
data_Z = [f reZ_cwt imZ_cwt];

%% taking imp data at specific frequency
for i = 1:74 %74
    time_plot(i) = time(sample_length*2*i)';
    %abs_plot(i) = abs(imp_cwt(end,sample_length*2*i));
    real_plot(i) = real(imp_cwt(end,sample_length*2*i)); %real_plot(i) = real(imp_cwt(end,sample_length*2*i)); real_plot(i) = real(imp_cwt(end-2,sample_length*i));
    %imag_plot(i) = imag(imp_cwt(1,sample_length*2*i));
end
time_plot_t = time_plot';
%abs_plot_t =  abs_plot';
real_plot_t = real_plot';
figure(14)
plot(time_plot, real_plot, '-ko');
%axis([0 75 2.2e5 3.3e5])
%% whole data set
resol = time(2)-time(1);
time_re_del = (1:1:length(real_del_cwt))'.*resol;
% figure(15)
% plot(time_re_del,real_imp_cwt(end-5,:),"k.")
%% creating data set
real_imp_cwt_t = real_imp_cwt';
%% creating data set for fitting
% f_t = transpose(f);
% imp_cwt_t = transpose(imp_cwt);
% for i = 1:74
%     Z_data(3*i-2:3*i,:) =  [f_t; real(imp_cwt_t(i*sample_length,:)); imag(imp_cwt_t(i*sample_length,:))]; %*sample_length
%     Z_data_t = transpose(Z_data);
% end
% %%
% options = optimoptions('fmincon','OptimalityTolerance',1e-11,'ConstraintTolerance',1e-11);
% 
% %% fitting
% for column_z = 1:74 %k
%     freq = Z_data_t(:,3*column_z-2);
%     reZ = Z_data_t(:,3*column_z-1);
%     imZ = Z_data_t(:,3*column_z);
%     % sample percentage and number of repetition
%  sample = 100;
%  N = 15;
%  sample_point=round(sample/100*length(f_t));
%  Z_raw=reZ + 1i*(-imZ);    
%     % assigning initial value of Randall circuit parameters
%  %max_imZ = max(imZ);
%  [Amax , idx] = max(imZ);
%  freq_max=freq(idx);
%  
%  Rs=  reZ(1);
%  Rct= reZ(length(f));
%  Cdl = 5e-10;
%  %Cdl = 1/(2*Rct*2*pi*freq_max);
%  n = 0.9;
% 
% initial=[Rs, Rct, Cdl,n]; 
% lower_bound=[0,0,0,0.7];
% upper_bound=[1e6, 5e7, 1e-8, 1];
% 
%     % getting sample points from EIS data
%     for index_N=1:N
%         index_point=randsample(length(f),sample_point);
%         index_point=sort(index_point);
%         sample_freq(:,index_N)=freq(index_point);
%         sample_Z(:,index_N)= Z_raw(index_point);
% 
%         % residue function
%         residue = @(x)sum(abs(sample_Z(:,index_N)-eq_circuit(x(1),x(2),x(3),x(4),sample_freq(:,index_N))).^2); %x(4),
%     
%         x(index_N,:)= fmincon(residue,initial,[],[],[],[],lower_bound,upper_bound,[],options); %[],options);
%     
%         fit_z(index_N,:)=eq_circuit(x(index_N,1), x(index_N,2), x(index_N,3),x(index_N,4), sample_freq(:,index_N)); % x(index_N,4),
%         
%         %fit_Z(:,column_z) = fit_z;
%     end   
%     
%     % fitting avarage
%     ave_Rs= mean(x(:,1));
%     ave_Rct= mean(x(:,2));
%     ave_Cdl= mean(x(:,3));
%     ave_n = mean(x(:,4));
%     fit_z_ave = eq_circuit(ave_Rs, ave_Rct, ave_Cdl, ave_n, freq);
%     fit_z_ave_1(:,column_z) = fit_z_ave(1:end);
%     ave_para = [ave_Rs ave_Rct ave_Cdl ave_n];
%     ave_para_1(:, column_z) = ave_para;
%     
%     % plot
% %     figure(30);
% %     plot(conj(fit_z_ave),'ko');
% %     hold on;
% %     axis equal;
% end
% hold off;
% ave_para_transpose = transpose(ave_para_1);   
% 
% %% compare fitting with raw
% 
% compare_no = 11;
% figure(22)
% plot(real(fit_z_ave_1(:,compare_no)) , imag(fit_z_ave_1(:,compare_no)), 'k-');
% hold on;
% plot(real(imp_cwt_t(compare_no*sample_length,:)), -imag(imp_cwt_t(compare_no*sample_length,:)), 'r.');
% hold off;
% %axis equal;
% %axis([0 inf 0 inf]);
% %%
% figure(23)
% for num_plot = 1:74 %k
%     plot3(arr+num_plot , real(fit_z_ave_1(:,num_plot)) , -imag(fit_z_ave_1(:,num_plot)) , 'k.-') ;
%     hold on;
%     grid on;
% end
% hold off;
%%  circuit function
function Z_eq=eq_circuit(rs,rct,cdl,n1,f) %W
    Xc=1./((1i*cdl*2*pi.*f).^n1);
    %Z_w=1./(W*sqrt(j*pi.*f));
    Z1=rct; %+Z_w;

    Z_eq=rs+Z1.*Xc./(Z1+Xc);
end