clear all;
path2=[ '../../DataSet/MonkeyECoG'];
addpath(genpath(path2));

N=128;

Tlag=4;
NLAG=8;

NUMWIN=20;

NPC=10;

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
        
        %% pca
        
        X=tss';
        X = bsxfun(@minus,X,mean(X));
        [coe,pcats,lat]=pca(X);
        tss2=pcats(:,1:NPC)';

        %%
        
        for i=1:NUMWIN
            ts=tss2(:,1+Tmax*(i-1):Tlag:Tmax+Tmax*(i-1));
            GCf=granger(ts,NLAG);
%             GCr=granger(flip(ts,2),NLAG);
%             Reference=((GCf(:)-GCr(:)).^2)';
            GCft=GCf';
            Reference=((GCf(:)-GCft(:)).^2)';
            FowRev{xx,nsub}=nanmean(Reference);
            nsub=nsub+1;
        end
    end
    NSUB(xx)=nsub-1;
end

for sub=1:NSUB(1)
    Wake(sub)=FowRev{1,sub};
end
for sub=1:NSUB(2)
    Sleep(sub)=FowRev{2,sub};
end
for sub=1:NSUB(3)
    KTW(sub)=FowRev{3,sub};
end
for sub=1:NSUB(4)
    KT(sub)=FowRev{4,sub};
end
for sub=1:NSUB(5)
    KTR(sub)=FowRev{5,sub};
end
for sub=1:NSUB(6)
    KTMDW(sub)=FowRev{6,sub};
end
for sub=1:NSUB(7)
    KTMD(sub)=FowRev{7,sub};
end
for sub=1:NSUB(8)
    KTMDR(sub)=FowRev{8,sub};
end
for sub=1:NSUB(9)
    PFW(sub)=FowRev{9,sub};
end
for sub=1:NSUB(10)
    PF(sub)=FowRev{10,sub};
end
for sub=1:NSUB(11)
    PFR(sub)=FowRev{11,sub};
end
for sub=1:NSUB(12)
    MDW(sub)=FowRev{12,sub};
end
for sub=1:NSUB(13)
    MD(sub)=FowRev{13,sub};
end
for sub=1:NSUB(14)
    MDR(sub)=FowRev{14,sub};
end

figure(1);
subplot(1,5,1)
violinplot([Wake Sleep],[zeros(1,NSUB(1)),ones(1,NSUB(2))]);
ranksum(Wake,Sleep)

subplot(1,5,5)
violinplot([KTW' KT' KTR']);
ranksum(KTW,KT)
ranksum(KT,KTR)
ranksum(KTW,KTR)

subplot(1,5,4)
violinplot([KTMDW KTMD KTMDR],[zeros(1,NSUB(6)),ones(1,NSUB(7)),2*ones(1,NSUB(8))]);
ranksum(KTMDW,KTMD)
ranksum(KTMD,KTMDR)
ranksum(KTMDW,KTMDR)

subplot(1,5,2)
violinplot([PFW' PF' PFR']);
ranksum(PFW,PF)
ranksum(PF,PFR)
ranksum(PFW,PFR)

subplot(1,5,3)
violinplot([MDW' MD' MDR']);
ranksum(MDW,MD)
ranksum(MD,MDR)
ranksum(MDW,MDR)


save results_GCchronosBlock_Monkey_asym.mat FowRev NSUB;
% save results_GCchronosBlock_Monkey.mat FowRev NSUB;

