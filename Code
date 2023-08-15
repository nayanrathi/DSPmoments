clc;
%clear all;
%pkg load statistics;
r_no=75;
C_Data=zeros(1,256);
M0_Data=zeros(1,r_no);
M1_Data=zeros(1,r_no);
SNR_Data=zeros(1,r_no);
SNR_Data1=zeros(1,r_no);
DopplerWidth=zeros(1,r_no);
count=0;
for xr=1:r_no
for range_bin=xr:xr;
  %range_bin=98;
    for k=5:8
               Framecount = 0;%This variable is used to store no of frames in the binary file
HEADERSIZE = 128; %Header size in Bytes in the binary file.
%--------------------------------------------------------------------------
% Data Reading
%--------------------------------------------------------------------------
%while(1)
%Get frame no to display from user
 %FrameNO =sscanf(input('Display Frame: ', 's'), '%d');
 Offset = 0;
 %Open binary file for reading
%fid(1) = fopen('F:\Matlab_Prgrms\matlab_c\13MA2017SHT1.r20','r');
fid(1) = fopen('C:\Users\sansk\Downloads\BE Project\24NO2021SHT1.r1','r');
%Read header
[A] = fread(fid(1),64,'int16');
%copy no of range bins and no of FFT poits
baudlength = A(2);
w1len = A(12);  %Window Length
nrgb = A(3);  % Number of Range Points
nfft = A(4);  % Number of FFT points
nci = A(5);  % Number of Coherent Integrations
nici = A(6);  % Number of Incoherent Integrations
ipp = A(7);  % Interpulse Period
pw = A(8);  % Pulse Width
w1start = A(11);  % Window Starts at at 10 micro seconds 1.5 km
nbeam = A(21);  % Number of Beams
%Loop to skip the headers up to user enterd frame header
 for n=1:(k-1)
    datasize = nrgb*nfft;
    Offset = Offset + (HEADERSIZE + (datasize*8)); %8: I/Q=2,Channels=1,Byte=4,
    fseek(fid(1),Offset, -1);
    [A] = fread(fid(1),64,'int16');
    nrgb = A(3);
    nfft = A(4);
 end


 for ft=1:nfft
    for rgb=1:nrgb
        I_Data_temp(rgb,ft) = fread(fid(1),1,'int32');
        Q_Data_temp(rgb,ft) = fread(fid(1),1,'int32');
    end
 end

%  Data_Frame=I_Data_temp+j*Q_Data_temp;
%  xlswrite('MST_Radar.xls',Data_Frame)

%Time domain data plotting
% Calculations of parameters for plotting
baud  = baudlength; %w1len/nrgb;                % Baud length
rrm   = ((3e8*(baud*1e-6))/2);    % Range resolution in mts
rrk   = rrm/1000;                 % Range resolution in kms
x     = 1/(ipp*1e-6*nci);        % Sampling frequnecy
res   = x/nfft;                   % Doppler resolution in Hz
xax   = linspace(-nfft/2,nfft/2-1,nfft); % X-axis FFT points
fn    = xax*res; % X-axis
wst   = (w1start)*rrk; % start height im kms (refer WPRS Timing diagram for  sampling window starts at starting of TX Pulse)
hh    = linspace(wst,(wst+(rrk*nrgb)), nrgb); % height in kms
td_x  = (ipp*nci*nfft)/1000;
%td_arranged = transpose(linspace(0,td_x,nfft));
td_arranged = (linspace(-3, 3,nfft));

%---------plotting of IQ Time domain data start here------------
for i = range_bin
    I_data = I_Data_temp(i,:);
    Q_data = Q_Data_temp(i,:);

    I_data=0.15*I_data/max(abs(I_data)); % this is to reduce the value for plotting
    Q_data=0.15*Q_data/max(abs(Q_data)); % this is to reduce the value for plotting

end

    end

 CMPLX_DATA=I_data+j*Q_data;

    Spectra=(abs(fft(CMPLX_DATA)));  %Fourier Spectrum for Incoherent Integration

    C_Data=C_Data+Spectra;
    end
    Power_Spectra=((circshift((C_Data/4),[0,128])).^2)./256;   %Power Spectrum
    Power_Spectra(129)=(Power_Spectra(128)+Power_Spectra(130))/2;  %Power Spectrum Cleaning
    %plot(Power_Spectra)

    %%Noise Level Estimation
Spectra_Sort=sort(Power_Spectra);
%%maxSn = nanmax(Power_Spectra);
%%N = length(Power_Spectra);   % Number of spectral estimates
%%Sn = sort(Power_Spectra, 'descend');
%%SnR2 = Sn(end); % Set signal-to-noise threshold to lowest value in spectrum.
%%for i = 1:N
%%    n = sum(isfinite(Sn(i:end)));
%%    P = nansum(Sn(i:end)) / n;          % Mean noise level for current iteration.
%%    Q = nansum(Sn(i:end).^2) / n - P^2; % Variance of noise level for current iteration.
%%    R2 = P^2 / (Q);
%%    if (R2 > 1)                         % White noise criteria has been met.
%%        SnR2 = Sn(i);
%%        break
%%    end
%%end


    for p=1:256
        P=0;
        Q=0;
        for k=1:p
        P=P+Spectra_Sort(k);
        Q=Q+((Spectra_Sort(k))^2);
        end
        Pn(p)=P/p;
        Qn(p)=Q/p;
        Qnn(p)=Qn(p)-(Pn(p))^2;
        R(p)=Pn(p)^2/(Qnn(p));
%         Qnn(p)=Qn(p)-(Pn(p))^2;
%         if Qnn(p)>0
%             R(p)=Pn(p)^2/(Qnn(p));
%         end
    end
 for l=1:256
     if R(l)>1
         K_s=l;
     end
 end

%%Calculation of Spectral Peak, Lower and Uppder Bounds, Doppler Velocity

Power_Spectra_L=Power_Spectra-Spectra_Sort(K_s);  % Noise Level Subtracted i,e Signal Spectra
%Power_Spectra_L=Power_Spectra-SnR2;  % Noise Level Subtracted i,e Signal Spectra

%plot(Power_Spectra_L)
[peak_Value , peak_Index] = max(Power_Spectra_L); %Dopper Peak Index
peak_Index_k=peak_Index;
Doppler_Frequency(range_bin)=((peak_Index_k-129)*0.0954);
Doppler_Velocity_peak=((peak_Index_k-129)*0.0954*3e8)/(2*53e6);
 KK=peak_Index;
     while Power_Spectra_L(KK)>=0   %Power_Spectra_L(KK-1)
         KK=KK-1;
     end
     K_min=KK;    %Doppler Lower Point from the Peak
     KK=peak_Index;
     while Power_Spectra_L(KK)>=0  %Power_Spectra_L(KK+1)
         KK=KK+1;
     end
     K_max=KK;     %Doppler Uppder Point from the Peak

     %Calculation of Moments
     Power_Mean=0;
     Moment_First=0;
     Moment_Second=0;
     Power_Mean1=0;

     for k=K_min:K_max
         Power_Mean=Power_Mean+Power_Spectra_L(k);  %Total Power in Doppler Spectrum
     end
%      for k=K_min:K_max
%          Power_Mean1=Power_Mean1+Power_Spectra(k);  %Total Power in Doppler Spectrum
%      end

     for k=K_min:K_max
         f_i(k)=((k-128)/((160e-6)*K_max*256));
         Moment_First=Moment_First+f_i(k)*Power_Spectra_L(k);  %Doppler Frequency
     end
     for k=K_min:K_max
         f_i(k)=((k-128)/((160e-6)*K_max*256));
         Moment_Second=Moment_Second+(Power_Spectra_L(k))*(f_i(k)-(Moment_First/Power_Mean))^2;
     end
      Doppler_Velocity(xr)=(1/Power_Mean)*((Moment_First/Power_Mean*3e8)/(2*53e6));
 count=count+1;
 disp(count);
M0_Data(xr)=Power_Mean;
M1_Data(xr)=Moment_First/Power_Mean;
SNR_Data(xr)=10*log(((Power_Mean)/(256*Spectra_Sort(K_s))));
DopplerWidth(xr)=sqrt(Moment_Second/Power_Mean);
Pss=0;
K_s1=K_s-1;
for k1=1:K_s1
Pss=Pss+Power_Spectra(k1);
end
Ps=Pss/(K_s1);
Pnn=0;
for k2=K_s:256
Pnn=Pnn+Power_Spectra(k2);
end
Pn=Pnn/(256-K_s);
SNR_Data1(xr)=Ps/Pn;
%%fprintf('Moment first: %d\n',Moment_First/Power_Mean);
%%fprintf('Moment Second: %d\n',Moment_Second);
%%plot(xr,Moment_First/Power_Mean);
%%hold on;
fclose('all');
end
% subplot(1,6,1);
% plot(M1_Data,(1:r_no));
% %xlim([-1 1]);
% title('Moments(mean doppler)');
% %xlim([-1 1]);
% subplot(1,6,2);
% plot(M0_Data,(1:r_no));
% title('Moments(signal power)');
%  %Doppler_Velocity=(1/Power_Mean)*(Moment_First*3e8)/(2*53e6);
% subplot(1,6,3);
% plot(SNR_Data,(1:r_no));
% title('SNR Ratio');
% subplot(1,6,4);
% plot(DopplerWidth,(1:r_no));
% title('Doppler Width');
% % subplot(1,6,5);
% % plot(SNR_Data1,(1:r_no));
% % title('SNR_Data Ps/Pn');
% subplot(1,6,6);
% plot(10*log(SNR_Data1),(1:r_no));
% title('SNR_Data Ps/Pn');
% Warning: Imaginary parts of complex X and/or Y arguments ignored 
subplot(1,5,1);
plot(M1_Data,(1:r_no));
%xlim([-1 1]);
title('Moments(mean doppler)');
%xlim([-1 1]);
subplot(1,5,2);
plot(M0_Data,(1:r_no));
title('Moments(signal power)');
 %Doppler_Velocity=(1/Power_Mean)*(Moment_First*3e8)/(2*53e6);
subplot(1,5,3);
plot(SNR_Data,(1:r_no));
title('SNR Ratio');
subplot(1,5,4);
plot(2*DopplerWidth,(1:r_no));
title('Doppler Width');
% subplot(1,6,5);
% plot(SNR_Data1,(1:r_no));
% title('SNR_Data Ps/Pn');
subplot(1,5,5);
plot(10*log(SNR_Data1),(1:r_no));
title('SNR_Data Ps/Pn');
%%     Doppler_Frequ_Dispersion=(1/Power_Mean)*Moment_Second;
%%     Doppler_Width=2*sqrt(Moment_Second);

% range_arranged=linspace(1,140,140)
%plot(Doppler_Velocity_peak(range_bin), range_arranged,'r*-')
%axis([-3 3 1.5 22.5]);


