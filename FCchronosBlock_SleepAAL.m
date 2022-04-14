clear all;
load DataSleepW_N3.mat;

%% Example for comparison of two conditions....

N=90;
NSUB=15;
Tau=1;

% Parameters of the data
TR=2.08;  % Repetition Time (seconds)
% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.008;                    % lowpass frequency of filter (Hz)
fhi = 0.08;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter


FowRev_W=zeros(1,NSUB);
FowRev_N3=zeros(1,NSUB);

for sub=1:NSUB  % over subjects
    sub
    ts=TS_W{sub};
    clear signal_filt;
    for seed=1:N
        ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));    %filtering
        signal_filt(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
    end
    ts=signal_filt(:,10:end-10);
    Tm=size(ts,2);
    FCtf=corr(ts(:,1:Tm-Tau)',ts(:,1+Tau:Tm)');       %% Core...FC tau foward
    FCtr=corr(ts(:,Tm:-1:Tau+1)',ts(:,Tm-Tau:-1:1)'); %% FC tau reversal
    Itauf=-0.5*log(1-FCtf.*FCtf);  %% Mutual information...
    Itaur=-0.5*log(1-FCtr.*FCtr);
    Reference=((Itauf(:)-Itaur(:)).^2)';
    index=find(Reference>quantile(Reference,0.0));
    FowRev_W(sub)=nanmean(Reference(index));
end

for sub=1:NSUB
    sub
    ts=TS_N3{sub};
    clear signal_filt;
    for seed=1:N
        ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
        signal_filt(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
    end
    ts=signal_filt(:,10:end-10);
    Tm=size(ts,2);
    FCtf=corr(ts(:,1:Tm-Tau)',ts(:,1+Tau:Tm)');
    FCtr=corr(ts(:,Tm:-1:Tau+1)',ts(:,Tm-Tau:-1:1)');
    Itauf=-0.5*log(1-FCtf.*FCtf);
    Itaur=-0.5*log(1-FCtr.*FCtr);
    Reference=((Itauf(:)-Itaur(:)).^2)';
    index=find(Reference>quantile(Reference,0.0));
    FowRev_N3(sub)=nanmean(Reference(index));
end

figure(1);
violinplot([FowRev_W' FowRev_N3']);
ranksum(FowRev_W,FowRev_N3)

save results_FCchronosBlock_SleepAAL.mat FowRev_W FowRev_N3;

