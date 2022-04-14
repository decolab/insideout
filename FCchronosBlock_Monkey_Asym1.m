clear all;
path2=[ '../../DataSet/MonkeyECoG'];
addpath(genpath(path2));

N=128;
Tau=4;

Chibi=0;
George=0;
Kin2=0;
Su=1;

if Chibi==1
    s1=[1 2 3 8];
    s2=[1 2 3 10 11]; 
    s3=[1 4];
    s4=[1 4];
    s5=[1 4];
    s6=[10 11];
    s7=[10 11];
    s8=[10 10];
    s9=[1 3];
    s10=[1 3];
    s11=[1 3];
    s12=[1 3];
    s13=[1 3];
    s14=[1 3];
end
if George==1
    s1=[4 5 6 7];
    s2=[4 5 6 7 8 9]; 
    s3=[2 3];
    s4=[2 3];
    s5=[2 3];
    s6=[1 2 3];
    s7=[1 2 3];
    s8=[1 2 3];
    s9=[2 4];
    s10=[2 4];
    s11=[2 4];
    s12=[2 4];
    s13=[2 4];
    s14=[2 4];
end
if Kin2==1
    s6=[4 6 7];
    s7=[4 6 7];
    s8=[4 6 7];
end
if Su==1
    s6=[5 8 9];
    s7=[5 8 9];
    s8=[5 8 9];
end

task={'Wake';'Sleep';'KTW';'KT';'KTR';'KTMDW';'KTMD';'KTMDR';'PFW';'PF';'PFR';'MDW';'MD';'MDR'};

for xx=1:size(task,1)
    xx
    list=dir(['../../DataSet/MonkeyECoG/' task{xx} '/*.mat']);
    NSUB(xx)=length(list);
    
    for sub=1:NSUB(xx)
        load(list(sub).name);
        ts=TimeSeries;
        clear signal_filt;
        for seed=1:N
            ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
        end
        Tm=size(ts,2);
        FCtf=corr(ts(:,1:Tm-Tau)',ts(:,1+Tau:Tm)');
        Itauf=-0.5*log(1-FCtf.*FCtf);
        Iout(sub,:)=mean(Itauf,2)';
        Iin(sub,:)=mean(Itauf);
    end
    if Chibi==1 || George==1
        if xx==1
            Wakef=squeeze(mean(Iout(s1,:)));
            Waker=squeeze(mean(Iin(s1,:)));
        end
        if xx==2
            Sleepf=squeeze(mean(Iout(s2,:)));
            Sleepr=squeeze(mean(Iin(s2,:)));
        end
        if xx==3
            KTWf=squeeze(mean(Iout(s3,:)));
            KTWr=squeeze(mean(Iin(s3,:)));
        end
        if xx==4
            KTf=squeeze(mean(Iout(s4,:)));
            KTr=squeeze(mean(Iin(s4,:)));
        end
        if xx==5
            KTRf=squeeze(mean(Iout(s5,:)));
            KTRr=squeeze(mean(Iin(s5,:)));
        end
    end
    if xx==6
        KTMDWf=squeeze(mean(Iout(s6,:)));
        KTMDWr=squeeze(mean(Iin(s6,:)));
    end
    if xx==7
        KTMDf=squeeze(mean(Iout(s7,:)));
        KTMDr=squeeze(mean(Iin(s7,:)));
    end
    if xx==8
        KTMDRf=squeeze(mean(Iout(s8,:)));
        KTMDRr=squeeze(mean(Iin(s8,:)));
    end
    if Chibi==1 || George==1
        if xx==9
            PFWf=squeeze(mean(Iout(s9,:)));
            PFWr=squeeze(mean(Iin(s9,:)));
        end
        if xx==10
            PFf=squeeze(mean(Iout(s10,:)));
            PFr=squeeze(mean(Iin(s10,:)));
        end
        if xx==11
            PFRf=squeeze(mean(Iout(s11,:)));
            PFRr=squeeze(mean(Iin(s11,:)));
        end
        if xx==12
            MDWf=squeeze(mean(Iout(s12,:)));
            MDWr=squeeze(mean(Iin(s12,:)));
        end
        if xx==13
            MDf=squeeze(mean(Iout(s13,:)));
            MDr=squeeze(mean(Iin(s13,:)));
        end
        if xx==14
            MDRf=squeeze(mean(Iout(s14,:)));
            MDRr=squeeze(mean(Iin(s14,:)));
        end
    end
end

if Chibi==1
    save results_asym_Chibi.mat Wakef Sleepf KTWf KTf KTRf KTMDWf KTMDf KTMDRf ...
        PFWf PFf PFRf MDWf MDf MDRf ...
        Waker Sleepr KTWr KTr KTRr KTMDWr KTMDr KTMDRr ...
        PFWr PFr PFRr MDWr MDr MDRr;
end
if George==1
    save results_asym_George.mat Wakef Sleepf KTWf KTf KTRf KTMDWf KTMDf KTMDRf ...
        PFWf PFf PFRf MDWf MDf MDRf ...
        Waker Sleepr KTWr KTr KTRr KTMDWr KTMDr KTMDRr ...
        PFWr PFr PFRr MDWr MDr MDRr;
end
if Kin2==1
    save results_asym_Kin2.mat KTMDWf KTMDf KTMDRf ...
        KTMDWr KTMDr KTMDRr;
end
if Su==1
    save results_asym_Su.mat KTMDWf KTMDf KTMDRf ...
        KTMDWr KTMDr KTMDRr;
end
