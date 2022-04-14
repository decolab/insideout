clear all;
path2=[ '../../DataSet/MonkeyECoG'];
addpath(genpath(path2));

N=128;
Tau=1;
NUMWIN=20;

task={'Wake';'Sleep';'KTW';'KT';'KTR';'KTMDW';'KTMD';'KTMDR';'PFW';'PF';'PFR';'MDW';'MD';'MDR'};

for xx=1:size(task,1)
    xx
    list=dir(['../../DataSet/MonkeyECoG/' task{xx} '/*.mat']);
    NSUB(xx)=length(list);
    
    nsub=1;
    for sub=1:NSUB(xx)
        sub
        load(list(sub).name);
        ts=TimeSeries;
        clear signal_filt;
        for seed=1:N
            ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
        end
        tss=ts;
        Tmm=size(tss,2);
        Tmax=floor(Tmm/NUMWIN);
        
        %%
        
        for i=1:NUMWIN
            ts=tss(:,1+Tmax*(i-1):Tmax+Tmax*(i-1));
            Tm=size(ts,2);
            FCtf=corr(ts(:,1:Tm-Tau)',ts(:,1+Tau:Tm)');
            Itauf=-0.5*log(1-FCtf.*FCtf);
            Itaufout=sum(Itauf);
            Itaufin=sum(Itauf,2)';
            for j=1:N
                HierNodeout{xx,nsub,j}=Itaufout(j);
                HierNodein{xx,nsub,j}=Itaufin(j);
            end
            nsub=nsub+1;
        end
    end
    NSUB(xx)=nsub-1;
end


for sub=1:NSUB(1)
    Wakeo(sub,:)=horzcat(HierNodeout{1,sub,:});
    Wakei(sub,:)=horzcat(HierNodein{1,sub,:});
end
Wakeout=squeeze(mean(Wakeo));
Wakein=squeeze(mean(Wakei));
Wakehier=0.5*(std(Wakeo)+std(Wakei));

for sub=1:NSUB(2)
    Sleepo(sub,:)=horzcat(HierNodeout{2,sub,:});
    Sleepi(sub,:)=horzcat(HierNodein{2,sub,:});
end
Sleepout=squeeze(mean(Sleepo));
Sleepin=squeeze(mean(Sleepi));
Sleephier=0.5*(std(Sleepo)+std(Sleepi));

for sub=1:NSUB(3)
    KTWo(sub,:)=horzcat(HierNodeout{3,sub,:});
    KTWi(sub,:)=horzcat(HierNodein{3,sub,:});
end
KTWout=squeeze(mean(KTWo));
KTWin=squeeze(mean(KTWi));
KTWhier=0.5*(std(KTWo)+std(KTWi));

for sub=1:NSUB(4)
    KTo(sub,:)=horzcat(HierNodeout{4,sub,:});
    KTi(sub,:)=horzcat(HierNodein{4,sub,:});
end
KTout=squeeze(mean(KTo));
KTin=squeeze(mean(KTi));
KThier=0.5*(std(KTo)+std(KTi));


for sub=1:NSUB(5)
    KTRo(sub,:)=horzcat(HierNodeout{5,sub,:});
    KTRi(sub,:)=horzcat(HierNodein{5,sub,:});
end
KTRout=squeeze(mean(KTRo));
KTRin=squeeze(mean(KTRi));
KTRhier=0.5*(std(KTRo)+std(KTRi));

for sub=1:NSUB(6)
    KTMDWo(sub,:)=horzcat(HierNodeout{6,sub,:});
    KTMDWi(sub,:)=horzcat(HierNodein{6,sub,:});
end
KTMDWout=squeeze(mean(KTMDWo));
KTMDWin=squeeze(mean(KTMDWi));
KTMDWhier=0.5*(std(KTMDWo)+std(KTMDWi));

for sub=1:NSUB(7)
    KTMDo(sub,:)=horzcat(HierNodeout{7,sub,:});
    KTMDi(sub,:)=horzcat(HierNodein{7,sub,:});
end
KTMDout=squeeze(mean(KTMDo));
KTMDin=squeeze(mean(KTMDi));
KTMDhier=0.5*(std(KTMDo)+std(KTMDi));


for sub=1:NSUB(8)
    KTMDRo(sub,:)=horzcat(HierNodeout{8,sub,:});
    KTMDRi(sub,:)=horzcat(HierNodein{8,sub,:});
end
KTMDRout=squeeze(mean(KTMDRo));
KTMDRin=squeeze(mean(KTMDRi));
KTMDRhier=0.5*(std(KTMDRo)+std(KTMDRi));

for sub=1:NSUB(9)
    PFWo(sub,:)=horzcat(HierNodeout{9,sub,:});
    PFWi(sub,:)=horzcat(HierNodein{9,sub,:});
end
PFWout=squeeze(mean(PFWo));
PFWin=squeeze(mean(PFWi));
PFWhier=0.5*(std(PFWo)+std(PFWi));

for sub=1:NSUB(10)
    PFo(sub,:)=horzcat(HierNodeout{10,sub,:});
    PFi(sub,:)=horzcat(HierNodein{10,sub,:});
end
PFout=squeeze(mean(PFo));
PFin=squeeze(mean(PFi));
PFhier=0.5*(std(PFo)+std(PFi));


for sub=1:NSUB(11)
    PFRo(sub,:)=horzcat(HierNodeout{11,sub,:});
    PFRi(sub,:)=horzcat(HierNodein{11,sub,:});
end
PFRout=squeeze(mean(PFRo));
PFRin=squeeze(mean(PFRi));
PFRhier=0.5*(std(PFRo)+std(PFRi));

for sub=1:NSUB(12)
    MDWo(sub,:)=horzcat(HierNodeout{12,sub,:});
    MDWi(sub,:)=horzcat(HierNodein{12,sub,:});
end
MDWout=squeeze(mean(MDWo));
MDWin=squeeze(mean(MDWi));
MDWhier=0.5*(std(MDWo)+std(MDWi));

for sub=1:NSUB(13)
    MDo(sub,:)=horzcat(HierNodeout{13,sub,:});
    MDi(sub,:)=horzcat(HierNodein{13,sub,:});
end
MDout=squeeze(mean(MDo));
MDin=squeeze(mean(MDi));
MDhier=0.5*(std(MDo)+std(MDi));

for sub=1:NSUB(14)
    MDRo(sub,:)=horzcat(HierNodeout{14,sub,:});
    MDRi(sub,:)=horzcat(HierNodein{14,sub,:});
end
MDRout=squeeze(mean(MDRo));
MDRin=squeeze(mean(MDRi));
MDRhier=0.5*(std(MDRo)+std(MDRi));

save results_FCchronosBlock_Monkey_FCtau.mat HierNodeout HierNodein NSUB ...
    Wakeout Wakein Wakehier Sleepout Sleepin Sleephier KTWout KTWin KTWhier KTout KTin KThier KTRout KTRin KTRhier ...
    KTMDWout KTMDWin KTMDWhier KTMDout KTMDin KTMDhier KTMDRout KTMDRin KTMDRhier ...
    PFWout PFWin PFWhier PFout PFin PFhier PFRout PFRin PFRhier ...
    MDWout MDWin MDWhier MDout MDin MDhier MDRout MDRin MDRhier;
 

