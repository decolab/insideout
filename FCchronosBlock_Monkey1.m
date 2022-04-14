clear all;
path2=[ '../../DataSet/MonkeyECoG'];
addpath(genpath(path2));

N=128;
Tau=4;
NUMWIN=20;

Chibi=1;
George=1;
Kin2=1;
Su=1;

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
        [coe(xx,sub,:,:),pcats,lat(xx,sub,:)]=pca(X);
        tss2=pcats(:,1:NPC)';
        
        %%
        
        for i=1:NUMWIN
            ts=tss2(:,1+Tmax*(i-1):Tmax+Tmax*(i-1));
            Tm=size(ts,2);
            FCtf=corr(ts(:,1:Tm-Tau)',ts(:,1+Tau:Tm)');
            FCtr=corr(ts(:,Tm:-1:Tau+1)',ts(:,Tm-Tau:-1:1)');
            Itauf=-0.5*log(1-FCtf.*FCtf);
            Itaur=-0.5*log(1-FCtr.*FCtr);
            Reference=((Itauf(:)-Itaur(:)).^2)';
            index=find(Reference>quantile(Reference,0.0));
            FowRev{xx,nsub}=nanmean(Reference(index));
            Idiff=(Itauf-Itaur).^2;
            Idiffnpca=0.5*(mean(Idiff)+mean(Idiff,2)');
            Ioutpca=mean(Itauf,2)';
            Iinpca=mean(Itauf);
            Hierarchy{xx,nsub}=nanstd(Reference(index)); %%std(Idiffnpca);
            for npc=1:NPC
                FowRevPCA{xx,nsub,npc}=Idiffnpca(npc);
                OutPCA{xx,nsub,npc}=Ioutpca(npc);
                InPCA{xx,nsub,npc}=Iinpca(npc);
            end
            %             Cxttauf=cov([ts(:,1:Tm-Tau)' ts(:,1+Tau:Tm)']);
            %             MIf=-0.5*logdet(Cxttauf);
            %             Cxttaur=cov([ts(:,Tm:-1:Tau+1)' ts(:,Tm-Tau:-1:1)']);
            %             MIr=-0.5*logdet(Cxttaur);
            %             FowRev{xx,nsub}=abs(MIf-MIr);
            nsub=nsub+1;
        end
    end
    NSUB(xx)=nsub-1;
end

THRLOW=5;
THRHIGH=90;

for sub=1:NSUB(1)
    Wake0(sub)=FowRev{1,sub};
    WakePCA(sub,:)=horzcat(FowRevPCA{1,sub,:});
    WakefPCA(sub,:)=horzcat(OutPCA{1,sub,:});
    WakerPCA(sub,:)=horzcat(InPCA{1,sub,:});
    Wakehier(sub)=Hierarchy{1,sub};
end
Wake=rmoutliers(Wake0,'percentiles',[THRLOW THRHIGH]);
NSUB2(1)=length(Wake);
for sub=1:NSUB(2)
    Sleep0(sub)=FowRev{2,sub};
    SleepPCA(sub,:)=horzcat(FowRevPCA{2,sub,:});
    SleepfPCA(sub,:)=horzcat(OutPCA{2,sub,:});
    SleeprPCA(sub,:)=horzcat(InPCA{2,sub,:});
    Sleephier(sub)=Hierarchy{2,sub};
end
Sleep=rmoutliers(Sleep0,'percentiles',[THRLOW THRHIGH]);
NSUB2(2)=length(Sleep);
for sub=1:NSUB(3)
    KTW0(sub)=FowRev{3,sub};
    KTWPCA(sub,:)=horzcat(FowRevPCA{3,sub,:});
    KTWfPCA(sub,:)=horzcat(OutPCA{3,sub,:});
    KTWrPCA(sub,:)=horzcat(InPCA{3,sub,:});
    KTWhier(sub)=Hierarchy{3,sub};
end
KTW=rmoutliers(KTW0,'percentiles',[THRLOW THRHIGH]);
NSUB2(3)=length(KTW);
for sub=1:NSUB(4)
    KT0(sub)=FowRev{4,sub};
    KTPCA(sub,:)=horzcat(FowRevPCA{4,sub,:});
    KTfPCA(sub,:)=horzcat(OutPCA{4,sub,:});
    KTrPCA(sub,:)=horzcat(InPCA{4,sub,:});
    KThier(sub)=Hierarchy{4,sub};
end
KT=rmoutliers(KT0,'percentiles',[THRLOW THRHIGH]);
NSUB2(4)=length(KT);
for sub=1:NSUB(5)
    KTR0(sub)=FowRev{5,sub};
    KTRPCA(sub,:)=horzcat(FowRevPCA{5,sub,:});
    KTRfPCA(sub,:)=horzcat(OutPCA{5,sub,:});
    KTRrPCA(sub,:)=horzcat(InPCA{5,sub,:});
    KTRhier(sub)=Hierarchy{5,sub};
end
KTR=rmoutliers(KTR0,'percentiles',[THRLOW THRHIGH]);
NSUB2(5)=length(KTR);
for sub=1:NSUB(6)
    KTMDW0(sub)=FowRev{6,sub};
    KTMDWPCA(sub,:)=horzcat(FowRevPCA{6,sub,:});
    KTMDWfPCA(sub,:)=horzcat(OutPCA{6,sub,:});
    KTMDWrPCA(sub,:)=horzcat(InPCA{6,sub,:});
    KTMDWhier(sub)=Hierarchy{6,sub};
end
KTMDW=rmoutliers(KTMDW0,'percentiles',[THRLOW THRHIGH]);
NSUB2(6)=length(KTMDW);
for sub=1:NSUB(7)
    KTMD0(sub)=FowRev{7,sub};
    KTMDPCA(sub,:)=horzcat(FowRevPCA{7,sub,:});
    KTMDfPCA(sub,:)=horzcat(OutPCA{7,sub,:});
    KTMDrPCA(sub,:)=horzcat(InPCA{7,sub,:});
    KTMDhier(sub)=Hierarchy{7,sub};
end
KTMD=rmoutliers(KTMD0,'percentiles',[THRLOW THRHIGH]);
NSUB2(7)=length(KTMD);
for sub=1:NSUB(8)
    KTMDR0(sub)=FowRev{8,sub};
    KTMDRPCA(sub,:)=horzcat(FowRevPCA{8,sub,:});
    KTMDRfPCA(sub,:)=horzcat(OutPCA{8,sub,:});
    KTMDRrPCA(sub,:)=horzcat(InPCA{8,sub,:});
    KTMDRhier(sub)=Hierarchy{8,sub};
end
KTMDR=rmoutliers(KTMDR0,'percentiles',[THRLOW THRHIGH]);
NSUB2(8)=length(KTMDR);
for sub=1:NSUB(9)
    PFW0(sub)=FowRev{9,sub};
    PFWPCA(sub,:)=horzcat(FowRevPCA{9,sub,:});
    PFWfPCA(sub,:)=horzcat(OutPCA{9,sub,:});
    PFWrPCA(sub,:)=horzcat(InPCA{9,sub,:});
    PFWhier(sub)=Hierarchy{9,sub};
end
PFW=rmoutliers(PFW0,'percentiles',[THRLOW THRHIGH]);
NSUB2(9)=length(PFW);
for sub=1:NSUB(10)
    PF0(sub)=FowRev{10,sub};
    PFPCA(sub,:)=horzcat(FowRevPCA{10,sub,:});
    PFfPCA(sub,:)=horzcat(OutPCA{10,sub,:});
    PFrPCA(sub,:)=horzcat(InPCA{10,sub,:});
    PFhier(sub)=Hierarchy{10,sub};
end
PF=rmoutliers(PF0,'percentiles',[THRLOW THRHIGH]);
NSUB2(10)=length(PF);
for sub=1:NSUB(11)
    PFR0(sub)=FowRev{11,sub};
    PFRPCA(sub,:)=horzcat(FowRevPCA{11,sub,:});
    PFRfPCA(sub,:)=horzcat(OutPCA{11,sub,:});
    PFRrPCA(sub,:)=horzcat(InPCA{11,sub,:});
    PFRhier(sub)=Hierarchy{11,sub};
end
PFR=rmoutliers(PFR0,'percentiles',[THRLOW THRHIGH]);
NSUB2(11)=length(PFR);
for sub=1:NSUB(12)
    MDW0(sub)=FowRev{12,sub};
    MDWPCA(sub,:)=horzcat(FowRevPCA{12,sub,:});
    MDWfPCA(sub,:)=horzcat(OutPCA{12,sub,:});
    MDWrPCA(sub,:)=horzcat(InPCA{12,sub,:});
    MDWhier(sub)=Hierarchy{12,sub};
end
MDW=rmoutliers(MDW0,'percentiles',[THRLOW THRHIGH]);
NSUB2(12)=length(MDW);
for sub=1:NSUB(13)
    MD0(sub)=FowRev{13,sub};
    MDPCA(sub,:)=horzcat(FowRevPCA{13,sub,:});
    MDfPCA(sub,:)=horzcat(OutPCA{13,sub,:});
    MDrPCA(sub,:)=horzcat(InPCA{13,sub,:});
    MDhier(sub)=Hierarchy{13,sub};
end
MD=rmoutliers(MD0,'percentiles',[THRLOW THRHIGH]);
NSUB2(13)=length(MD);
for sub=1:NSUB(14)
    MDR0(sub)=FowRev{14,sub};
    MDRPCA(sub,:)=horzcat(FowRevPCA{14,sub,:});
    MDRfPCA(sub,:)=horzcat(OutPCA{14,sub,:});
    MDRrPCA(sub,:)=horzcat(InPCA{14,sub,:});
    MDRhier(sub)=Hierarchy{14,sub};
end
MDR=rmoutliers(MDR0,'percentiles',[THRLOW THRHIGH]);
NSUB2(14)=length(MDR);

figure(1);
violinplot([Wake Sleep],[zeros(1,NSUB2(1)),ones(1,NSUB2(2))]);
ranksum(Wake,Sleep)

figure(2);
violinplot([KTW' KT' KTR']);
ranksum(KTW,KT)
ranksum(KT,KTR)
ranksum(KTW,KTR)

figure(3);
violinplot([KTMDW KTMD KTMDR],[zeros(1,NSUB2(6)),ones(1,NSUB2(7)),2*ones(1,NSUB2(8))]);
ranksum(KTMDW,KTMD)
ranksum(KTMD,KTMDR)
ranksum(KTMDW,KTMDR)

figure(4);
violinplot([PFW' PF' PFR']);
ranksum(PFW,PF)
ranksum(PF,PFR)
ranksum(PFW,PFR)

figure(5);
violinplot([MDW' MD' MDR']);
ranksum(MDW,MD)
ranksum(MD,MDR)
ranksum(MDW,MDR)

%% Rendering Chibi
if Chibi==1
    clear WakeA SleepA KTWA KTA KTRA KTMDWA KTMDA KTMDRA ...
    PFWA PFA PFRA MDWA MDA MDRA ...
    WakefA SleepfA KTWfA KTfA KTRfA KTMDWfA KTMDfA KTMDRfA ...
    PFWfA PFfA PFRfA MDWfA MDfA MDRfA ...
    WakerA SleeprA KTWrA KTrA KTRrA KTMDWrA KTMDrA KTMDRrA ...
    PFWrA PFrA PFRrA MDWrA MDrA MDRrA;
    
    n=1;
    for s=[1 2 3 8]
        pc=squeeze(coe(1,s,:,:));
        we=mean(WakePCA(1+20*(s-1):20+20*(s-1),:));
        wef=mean(WakefPCA(1+20*(s-1):20+20*(s-1),:));
        wer=mean(WakerPCA(1+20*(s-1):20+20*(s-1),:));
        for i=1:NPC
            WakeA(n,i,:)=pc(:,i)*we(i);
            WakefA(n,i,:)=pc(:,i)*wef(i);
            WakerA(n,i,:)=pc(:,i)*wer(i);
        end
        n=n+1;
    end
    WakeA=mean(squeeze(mean(WakeA)));
    WakefA=mean(squeeze(mean(WakefA)));
    WakerA=mean(squeeze(mean(WakerA)));
    
    n=1;
    for s=[1 2 3 10 11]
        pc=squeeze(coe(2,s,:,:));
        we=mean(SleepPCA(1+20*(s-1):20+20*(s-1),:));
        wef=mean(SleepfPCA(1+20*(s-1):20+20*(s-1),:));
        wer=mean(SleeprPCA(1+20*(s-1):20+20*(s-1),:));
        for i=1:NPC
            SleepA(n,i,:)=pc(:,i)*we(i);
            SleepfA(n,i,:)=pc(:,i)*wef(i);
            SleeprA(n,i,:)=pc(:,i)*wer(i);
        end
        n=n+1;
    end
    SleepA=mean(squeeze(mean(SleepA)));
    SleepfA=mean(squeeze(mean(SleepfA)));
    SleeprA=mean(squeeze(mean(SleeprA)));
    
    n=1;
    for s=[1 4]
        pc=squeeze(coe(3,s,:,:));
        we=mean(KTWPCA(1+20*(s-1):20+20*(s-1),:));
        wef=mean(KTWfPCA(1+20*(s-1):20+20*(s-1),:));
        wer=mean(KTWrPCA(1+20*(s-1):20+20*(s-1),:));
        for i=1:NPC
            KTWA(n,i,:)=pc(:,i)*we(i);
            KTWfA(n,i,:)=pc(:,i)*wef(i);
            KTWrA(n,i,:)=pc(:,i)*wer(i);
        end
        n=n+1;
    end
    KTWA=mean(squeeze(mean(KTWA)));
    KTWfA=mean(squeeze(mean(KTWfA)));
    KTWrA=mean(squeeze(mean(KTWrA)));
    
    n=1;
    for s=[1 4]
        pc=squeeze(coe(4,s,:,:));
        we=mean(KTPCA(1+20*(s-1):20+20*(s-1),:));
        wef=mean(KTfPCA(1+20*(s-1):20+20*(s-1),:));
        wer=mean(KTrPCA(1+20*(s-1):20+20*(s-1),:));
        for i=1:NPC
            KTA(n,i,:)=pc(:,i)*we(i);
            KTfA(n,i,:)=pc(:,i)*wef(i);
            KTrA(n,i,:)=pc(:,i)*wer(i);
        end
        n=n+1;
    end
    KTA=mean(squeeze(mean(KTA)));
    KTfA=mean(squeeze(mean(KTfA)));
    KTrA=mean(squeeze(mean(KTrA)));
    
    n=1;
    for s=[1 4]
        pc=squeeze(coe(5,s,:,:));
        we=mean(KTRPCA(1+20*(s-1):20+20*(s-1),:));
        wef=mean(KTRfPCA(1+20*(s-1):20+20*(s-1),:));
        wer=mean(KTRrPCA(1+20*(s-1):20+20*(s-1),:));
        for i=1:NPC
            KTRA(n,i,:)=pc(:,i)*we(i);
            KTRfA(n,i,:)=pc(:,i)*wef(i);
            KTRrA(n,i,:)=pc(:,i)*wer(i);
        end
        n=n+1;
    end
    KTRA=mean(squeeze(mean(KTRA)));
    KTRfA=mean(squeeze(mean(KTRfA)));
    KTRrA=mean(squeeze(mean(KTRrA)));
    
    n=1;
    for s=[10 11]
        pc=squeeze(coe(6,s,:,:));
        we=mean(KTMDWPCA(1+20*(s-1):20+20*(s-1),:));
        wef=mean(KTMDWfPCA(1+20*(s-1):20+20*(s-1),:));
        wer=mean(KTMDWrPCA(1+20*(s-1):20+20*(s-1),:));
        for i=1:NPC
            KTMDWA(n,i,:)=pc(:,i)*we(i);
            KTMDWfA(n,i,:)=pc(:,i)*wef(i);
            KTMDWrA(n,i,:)=pc(:,i)*wer(i);
        end
        n=n+1;
    end
    KTMDWA=mean(squeeze(mean(KTMDWA)));
    KTMDWfA=mean(squeeze(mean(KTMDWfA)));
    KTMDWrA=mean(squeeze(mean(KTMDWrA)));
    
    n=1;
    for s=[10 11]
        pc=squeeze(coe(7,s,:,:));
        we=mean(KTMDPCA(1+20*(s-1):20+20*(s-1),:));
        wef=mean(KTMDfPCA(1+20*(s-1):20+20*(s-1),:));
        wer=mean(KTMDrPCA(1+20*(s-1):20+20*(s-1),:));
        for i=1:NPC
            KTMDA(n,i,:)=pc(:,i)*we(i);
            KTMDfA(n,i,:)=pc(:,i)*wef(i);
            KTMDrA(n,i,:)=pc(:,i)*wer(i);
        end
        n=n+1;
    end
    KTMDA=mean(squeeze(mean(KTMDA)));
    KTMDfA=mean(squeeze(mean(KTMDfA)));
    KTMDrA=mean(squeeze(mean(KTMDrA)));
    
    n=1;
    for s=[10 10]
        pc=squeeze(coe(8,s,:,:));
        we=mean(KTMDRPCA(1+20*(s-1):20+20*(s-1),:));
        wef=mean(KTMDRfPCA(1+20*(s-1):20+20*(s-1),:));
        wer=mean(KTMDRrPCA(1+20*(s-1):20+20*(s-1),:));
        for i=1:NPC
            KTMDRA(n,i,:)=pc(:,i)*we(i);
            KTMDRfA(n,i,:)=pc(:,i)*wef(i);
            KTMDRrA(n,i,:)=pc(:,i)*wer(i);
        end
        n=n+1;
    end
    KTMDRA=mean(squeeze(mean(KTMDRA)));
    KTMDRfA=mean(squeeze(mean(KTMDRfA)));
    KTMDRrA=mean(squeeze(mean(KTMDRrA)));
    
    n=1;
    for s=[1 3]
        pc=squeeze(coe(9,s,:,:));
        we=mean(PFWPCA(1+20*(s-1):20+20*(s-1),:));
        wef=mean(PFWfPCA(1+20*(s-1):20+20*(s-1),:));
        wer=mean(PFWrPCA(1+20*(s-1):20+20*(s-1),:));
        for i=1:NPC
            PFWA(n,i,:)=pc(:,i)*we(i);
            PFWfA(n,i,:)=pc(:,i)*wef(i);
            PFWrA(n,i,:)=pc(:,i)*wer(i);
        end
        n=n+1;
    end
    PFWA=mean(squeeze(mean(PFWA)));
    PFWfA=mean(squeeze(mean(PFWfA)));
    PFWrA=mean(squeeze(mean(PFWrA)));
    
    n=1;
    for s=[1 3]
        pc=squeeze(coe(10,s,:,:));
        we=mean(PFPCA(1+20*(s-1):20+20*(s-1),:));
        wef=mean(PFfPCA(1+20*(s-1):20+20*(s-1),:));
        wer=mean(PFrPCA(1+20*(s-1):20+20*(s-1),:));
        for i=1:NPC
            PFA(n,i,:)=pc(:,i)*we(i);
            PFfA(n,i,:)=pc(:,i)*wef(i);
            PFrA(n,i,:)=pc(:,i)*wer(i);
        end
        n=n+1;
    end
    PFA=mean(squeeze(mean(PFA)));
    PFfA=mean(squeeze(mean(PFfA)));
    PFrA=mean(squeeze(mean(PFrA)));
    
    n=1;
    for s=[1 3]
        pc=squeeze(coe(11,s,:,:));
        we=mean(PFRPCA(1+20*(s-1):20+20*(s-1),:));
        wef=mean(PFRfPCA(1+20*(s-1):20+20*(s-1),:));
        wer=mean(PFRrPCA(1+20*(s-1):20+20*(s-1),:));
        for i=1:NPC
            PFRA(n,i,:)=pc(:,i)*we(i);
            PFRfA(n,i,:)=pc(:,i)*wef(i);
            PFRrA(n,i,:)=pc(:,i)*wer(i);
        end
        n=n+1;
    end
    PFRA=mean(squeeze(mean(PFRA)));
    PFRfA=mean(squeeze(mean(PFRfA)));
    PFRrA=mean(squeeze(mean(PFRrA)));
    
    n=1;
    for s=[1 3]
        pc=squeeze(coe(12,s,:,:));
        we=mean(MDWPCA(1+20*(s-1):20+20*(s-1),:));
        wef=mean(MDWfPCA(1+20*(s-1):20+20*(s-1),:));
        wer=mean(MDWrPCA(1+20*(s-1):20+20*(s-1),:));
        for i=1:NPC
            MDWA(n,i,:)=pc(:,i)*we(i);
            MDWfA(n,i,:)=pc(:,i)*wef(i);
            MDWrA(n,i,:)=pc(:,i)*wer(i);
        end
        n=n+1;
    end
    MDWA=mean(squeeze(mean(MDWA)));
    MDWfA=mean(squeeze(mean(MDWfA)));
    MDWrA=mean(squeeze(mean(MDWrA)));
    
    n=1;
    for s=[1 3]
        pc=squeeze(coe(13,s,:,:));
        we=mean(MDPCA(1+20*(s-1):20+20*(s-1),:));
        wef=mean(MDfPCA(1+20*(s-1):20+20*(s-1),:));
        wer=mean(MDrPCA(1+20*(s-1):20+20*(s-1),:));
        for i=1:NPC
            MDA(n,i,:)=pc(:,i)*we(i);
            MDfA(n,i,:)=pc(:,i)*wef(i);
            MDrA(n,i,:)=pc(:,i)*wer(i);
        end
        n=n+1;
    end
    MDA=mean(squeeze(mean(MDA)));
    MDfA=mean(squeeze(mean(MDfA)));
    MDrA=mean(squeeze(mean(MDrA)));
    
    n=1;
    for s=[1 3]
        pc=squeeze(coe(14,s,:,:));
        we=mean(MDRPCA(1+20*(s-1):20+20*(s-1),:));
        wef=mean(MDRfPCA(1+20*(s-1):20+20*(s-1),:));
        wer=mean(MDRrPCA(1+20*(s-1):20+20*(s-1),:));
        for i=1:NPC
            MDRA(n,i,:)=pc(:,i)*we(i);
            MDRfA(n,i,:)=pc(:,i)*wef(i);
            MDRrA(n,i,:)=pc(:,i)*wer(i);
        end
        n=n+1;
    end
    MDRA=mean(squeeze(mean(MDRA)));
    MDRfA=mean(squeeze(mean(MDRfA)));
    MDRrA=mean(squeeze(mean(MDRrA)));
    
    save results_FCchronosBlock_Monkey_Chibi.mat FowRev FowRevPCA NSUB coe lat ...
    WakePCA SleepPCA KTWPCA KTPCA KTRPCA KTMDWPCA KTMDPCA KTMDRPCA ...
    PFWPCA PFPCA PFRPCA MDWPCA MDPCA MDRPCA ...
    WakeA SleepA KTWA KTA KTRA KTMDWA KTMDA KTMDRA ...
    PFWA PFA PFRA MDWA MDA MDRA ...
    WakefA SleepfA KTWfA KTfA KTRfA KTMDWfA KTMDfA KTMDRfA ...
    PFWfA PFfA PFRfA MDWfA MDfA MDRfA ...
    WakerA SleeprA KTWrA KTrA KTRrA KTMDWrA KTMDrA KTMDRrA ...
    PFWrA PFrA PFRrA MDWrA MDrA MDRrA ...
    Wakehier Sleephier KTWhier KThier KTRhier KTMDWhier KTMDhier KTMDRhier ...
    PFWhier PFhier PFRhier MDWhier MDhier MDRhier;
    
end
%%

%% Rendering George

if George==1
    clear WakeA SleepA KTWA KTA KTRA KTMDWA KTMDA KTMDRA ...
    PFWA PFA PFRA MDWA MDA MDRA ...
    WakefA SleepfA KTWfA KTfA KTRfA KTMDWfA KTMDfA KTMDRfA ...
    PFWfA PFfA PFRfA MDWfA MDfA MDRfA ...
    WakerA SleeprA KTWrA KTrA KTRrA KTMDWrA KTMDrA KTMDRrA ...
    PFWrA PFrA PFRrA MDWrA MDrA MDRrA;
    
    n=1;
    for s=[4 5 6 7]
        pc=squeeze(coe(1,s,:,:));
        we=mean(WakePCA(1+20*(s-1):20+20*(s-1),:));
        wef=mean(WakefPCA(1+20*(s-1):20+20*(s-1),:));
        wer=mean(WakerPCA(1+20*(s-1):20+20*(s-1),:));
        for i=1:NPC
            WakeA(n,i,:)=pc(:,i)*we(i);
            WakefA(n,i,:)=pc(:,i)*wef(i);
            WakerA(n,i,:)=pc(:,i)*wer(i);
        end
        n=n+1;
    end
    WakeA=mean(squeeze(mean(WakeA)));
    WakefA=mean(squeeze(mean(WakefA)));
    WakerA=mean(squeeze(mean(WakerA)));
    
    n=1;
    for s=[4 5 6 7 8 9]
        pc=squeeze(coe(2,s,:,:));
        we=mean(SleepPCA(1+20*(s-1):20+20*(s-1),:));
        wef=mean(SleepfPCA(1+20*(s-1):20+20*(s-1),:));
        wer=mean(SleeprPCA(1+20*(s-1):20+20*(s-1),:));
        for i=1:NPC
            SleepA(n,i,:)=pc(:,i)*we(i);
            SleepfA(n,i,:)=pc(:,i)*wef(i);
            SleeprA(n,i,:)=pc(:,i)*wer(i);
        end
        n=n+1;
    end
    SleepA=mean(squeeze(mean(SleepA)));
    SleepfA=mean(squeeze(mean(SleepfA)));
    SleeprA=mean(squeeze(mean(SleeprA)));
    
    n=1;
    for s=[2 3]
        pc=squeeze(coe(3,s,:,:));
        we=mean(KTWPCA(1+20*(s-1):20+20*(s-1),:));
        wef=mean(KTWfPCA(1+20*(s-1):20+20*(s-1),:));
        wer=mean(KTWrPCA(1+20*(s-1):20+20*(s-1),:));
        for i=1:NPC
            KTWA(n,i,:)=pc(:,i)*we(i);
            KTWfA(n,i,:)=pc(:,i)*wef(i);
            KTWrA(n,i,:)=pc(:,i)*wer(i);
        end
        n=n+1;
    end
    KTWA=mean(squeeze(mean(KTWA)));
    KTWfA=mean(squeeze(mean(KTWfA)));
    KTWrA=mean(squeeze(mean(KTWrA)));
    
    n=1;
    for s=[2 3]
        pc=squeeze(coe(4,s,:,:));
        we=mean(KTPCA(1+20*(s-1):20+20*(s-1),:));
        wef=mean(KTfPCA(1+20*(s-1):20+20*(s-1),:));
        wer=mean(KTrPCA(1+20*(s-1):20+20*(s-1),:));
        for i=1:NPC
            KTA(n,i,:)=pc(:,i)*we(i);
            KTfA(n,i,:)=pc(:,i)*wef(i);
            KTrA(n,i,:)=pc(:,i)*wer(i);
        end
        n=n+1;
    end
    KTA=mean(squeeze(mean(KTA)));
    KTfA=mean(squeeze(mean(KTfA)));
    KTrA=mean(squeeze(mean(KTrA)));
    
    n=1;
    for s=[2 3]
        pc=squeeze(coe(5,s,:,:));
        we=mean(KTRPCA(1+20*(s-1):20+20*(s-1),:));
        wef=mean(KTRfPCA(1+20*(s-1):20+20*(s-1),:));
        wer=mean(KTRrPCA(1+20*(s-1):20+20*(s-1),:));
        for i=1:NPC
            KTRA(n,i,:)=pc(:,i)*we(i);
            KTRfA(n,i,:)=pc(:,i)*wef(i);
            KTRrA(n,i,:)=pc(:,i)*wer(i);
        end
        n=n+1;
    end
    KTRA=mean(squeeze(mean(KTRA)));
    KTRfA=mean(squeeze(mean(KTRfA)));
    KTRrA=mean(squeeze(mean(KTRrA)));
    
    n=1;
    for s=[1 2 3]
        pc=squeeze(coe(6,s,:,:));
        we=mean(KTMDWPCA(1+20*(s-1):20+20*(s-1),:));
        wef=mean(KTMDWfPCA(1+20*(s-1):20+20*(s-1),:));
        wer=mean(KTMDWrPCA(1+20*(s-1):20+20*(s-1),:));
        for i=1:NPC
            KTMDWA(n,i,:)=pc(:,i)*we(i);
            KTMDWfA(n,i,:)=pc(:,i)*wef(i);
            KTMDWrA(n,i,:)=pc(:,i)*wer(i);
        end
        n=n+1;
    end
    KTMDWA=mean(squeeze(mean(KTMDWA)));
    KTMDWfA=mean(squeeze(mean(KTMDWfA)));
    KTMDWrA=mean(squeeze(mean(KTMDWrA)));
    
    n=1;
    for s=[1 2 3]
        pc=squeeze(coe(7,s,:,:));
        we=mean(KTMDPCA(1+20*(s-1):20+20*(s-1),:));
        wef=mean(KTMDfPCA(1+20*(s-1):20+20*(s-1),:));
        wer=mean(KTMDrPCA(1+20*(s-1):20+20*(s-1),:));
        for i=1:NPC
            KTMDA(n,i,:)=pc(:,i)*we(i);
            KTMDfA(n,i,:)=pc(:,i)*wef(i);
            KTMDrA(n,i,:)=pc(:,i)*wer(i);
        end
        n=n+1;
    end
    KTMDA=mean(squeeze(mean(KTMDA)));
    KTMDfA=mean(squeeze(mean(KTMDfA)));
    KTMDrA=mean(squeeze(mean(KTMDrA)));
    
    n=1;
    for s=[1 2 3]
        pc=squeeze(coe(8,s,:,:));
        we=mean(KTMDRPCA(1+20*(s-1):20+20*(s-1),:));
        wef=mean(KTMDRfPCA(1+20*(s-1):20+20*(s-1),:));
        wer=mean(KTMDRrPCA(1+20*(s-1):20+20*(s-1),:));
        for i=1:NPC
            KTMDRA(n,i,:)=pc(:,i)*we(i);
            KTMDRfA(n,i,:)=pc(:,i)*wef(i);
            KTMDRrA(n,i,:)=pc(:,i)*wer(i);
        end
        n=n+1;
    end
    KTMDRA=mean(squeeze(mean(KTMDRA)));
    KTMDRfA=mean(squeeze(mean(KTMDRfA)));
    KTMDRrA=mean(squeeze(mean(KTMDRrA)));
    
    n=1;
    for s=[2 4]
        pc=squeeze(coe(9,s,:,:));
        we=mean(PFWPCA(1+20*(s-1):20+20*(s-1),:));
        wef=mean(PFWfPCA(1+20*(s-1):20+20*(s-1),:));
        wer=mean(PFWrPCA(1+20*(s-1):20+20*(s-1),:));
        for i=1:NPC
            PFWA(n,i,:)=pc(:,i)*we(i);
            PFWfA(n,i,:)=pc(:,i)*wef(i);
            PFWrA(n,i,:)=pc(:,i)*wer(i);
        end
        n=n+1;
    end
    PFWA=mean(squeeze(mean(PFWA)));
    PFWfA=mean(squeeze(mean(PFWfA)));
    PFWrA=mean(squeeze(mean(PFWrA)));
    
    n=1;
    for s=[2 4]
        pc=squeeze(coe(10,s,:,:));
        we=mean(PFPCA(1+20*(s-1):20+20*(s-1),:));
        wef=mean(PFfPCA(1+20*(s-1):20+20*(s-1),:));
        wer=mean(PFrPCA(1+20*(s-1):20+20*(s-1),:));
        for i=1:NPC
            PFA(n,i,:)=pc(:,i)*we(i);
            PFfA(n,i,:)=pc(:,i)*wef(i);
            PFrA(n,i,:)=pc(:,i)*wer(i);
        end
        n=n+1;
    end
    PFA=mean(squeeze(mean(PFA)));
    PFfA=mean(squeeze(mean(PFfA)));
    PFrA=mean(squeeze(mean(PFrA)));
    
    n=1;
    for s=[2 4]
        pc=squeeze(coe(11,s,:,:));
        we=mean(PFRPCA(1+20*(s-1):20+20*(s-1),:));
        wef=mean(PFRfPCA(1+20*(s-1):20+20*(s-1),:));
        wer=mean(PFRrPCA(1+20*(s-1):20+20*(s-1),:));
        for i=1:NPC
            PFRA(n,i,:)=pc(:,i)*we(i);
            PFRfA(n,i,:)=pc(:,i)*wef(i);
            PFRrA(n,i,:)=pc(:,i)*wer(i);
        end
        n=n+1;
    end
    PFRA=mean(squeeze(mean(PFRA)));
    PFRfA=mean(squeeze(mean(PFRfA)));
    PFRrA=mean(squeeze(mean(PFRrA)));
    
    n=1;
    for s=[2 4]
        pc=squeeze(coe(12,s,:,:));
        we=mean(MDWPCA(1+20*(s-1):20+20*(s-1),:));
        wef=mean(MDWfPCA(1+20*(s-1):20+20*(s-1),:));
        wer=mean(MDWrPCA(1+20*(s-1):20+20*(s-1),:));
        for i=1:NPC
            MDWA(n,i,:)=pc(:,i)*we(i);
            MDWfA(n,i,:)=pc(:,i)*wef(i);
            MDWrA(n,i,:)=pc(:,i)*wer(i);
        end
        n=n+1;
    end
    MDWA=mean(squeeze(mean(MDWA)));
    MDWfA=mean(squeeze(mean(MDWfA)));
    MDWrA=mean(squeeze(mean(MDWrA)));
    
    n=1;
    for s=[2 4]
        pc=squeeze(coe(13,s,:,:));
        we=mean(MDPCA(1+20*(s-1):20+20*(s-1),:));
        wef=mean(MDfPCA(1+20*(s-1):20+20*(s-1),:));
        wer=mean(MDrPCA(1+20*(s-1):20+20*(s-1),:));
        for i=1:NPC
            MDA(n,i,:)=pc(:,i)*we(i);
            MDfA(n,i,:)=pc(:,i)*wef(i);
            MDrA(n,i,:)=pc(:,i)*wer(i);
        end
        n=n+1;
    end
    MDA=mean(squeeze(mean(MDA)));
    MDfA=mean(squeeze(mean(MDfA)));
    MDrA=mean(squeeze(mean(MDrA)));
    
    n=1;
    for s=[2 4]
        pc=squeeze(coe(14,s,:,:));
        we=mean(MDRPCA(1+20*(s-1):20+20*(s-1),:));
        wef=mean(MDRfPCA(1+20*(s-1):20+20*(s-1),:));
        wer=mean(MDRrPCA(1+20*(s-1):20+20*(s-1),:));
        for i=1:NPC
            MDRA(n,i,:)=pc(:,i)*we(i);
            MDRfA(n,i,:)=pc(:,i)*wef(i);
            MDRrA(n,i,:)=pc(:,i)*wer(i);
        end
        n=n+1;
    end
    MDRA=mean(squeeze(mean(MDRA)));
    MDRfA=mean(squeeze(mean(MDRfA)));
    MDRrA=mean(squeeze(mean(MDRrA)));
    
    save results_FCchronosBlock_Monkey_George.mat FowRev FowRevPCA NSUB coe lat ...
    WakePCA SleepPCA KTWPCA KTPCA KTRPCA KTMDWPCA KTMDPCA KTMDRPCA ...
    PFWPCA PFPCA PFRPCA MDWPCA MDPCA MDRPCA ...
    WakeA SleepA KTWA KTA KTRA KTMDWA KTMDA KTMDRA ...
    PFWA PFA PFRA MDWA MDA MDRA ...
    WakefA SleepfA KTWfA KTfA KTRfA KTMDWfA KTMDfA KTMDRfA ...
    PFWfA PFfA PFRfA MDWfA MDfA MDRfA ...
    WakerA SleeprA KTWrA KTrA KTRrA KTMDWrA KTMDrA KTMDRrA ...
    PFWrA PFrA PFRrA MDWrA MDrA MDRrA ...
    Wakehier Sleephier KTWhier KThier KTRhier KTMDWhier KTMDhier KTMDRhier ...
    PFWhier PFhier PFRhier MDWhier MDhier MDRhier;
    
end
%%

%% Rendering Kin2

if Kin2==1
    clear WakeA SleepA KTWA KTA KTRA KTMDWA KTMDA KTMDRA ...
    PFWA PFA PFRA MDWA MDA MDRA ...
    WakefA SleepfA KTWfA KTfA KTRfA KTMDWfA KTMDfA KTMDRfA ...
    PFWfA PFfA PFRfA MDWfA MDfA MDRfA ...
    WakerA SleeprA KTWrA KTrA KTRrA KTMDWrA KTMDrA KTMDRrA ...
    PFWrA PFrA PFRrA MDWrA MDrA MDRrA;   
    n=1;
    for s=[4 6 7]
        pc=squeeze(coe(6,s,:,:));
        we=mean(KTMDWPCA(1+20*(s-1):20+20*(s-1),:));
        wef=mean(KTMDWfPCA(1+20*(s-1):20+20*(s-1),:));
        wer=mean(KTMDWrPCA(1+20*(s-1):20+20*(s-1),:));
        for i=1:NPC
            KTMDWA(n,i,:)=pc(:,i)*we(i);
            KTMDWfA(n,i,:)=pc(:,i)*wef(i);
            KTMDWrA(n,i,:)=pc(:,i)*wer(i);
        end
        n=n+1;
    end
    KTMDWA=mean(squeeze(mean(KTMDWA)));
    KTMDWfA=mean(squeeze(mean(KTMDWfA)));
    KTMDWrA=mean(squeeze(mean(KTMDWrA)));
    
    n=1;
    for s=[4 6 7]
        pc=squeeze(coe(7,s,:,:));
        we=mean(KTMDPCA(1+20*(s-1):20+20*(s-1),:));
        wef=mean(KTMDfPCA(1+20*(s-1):20+20*(s-1),:));
        wer=mean(KTMDrPCA(1+20*(s-1):20+20*(s-1),:));
        for i=1:NPC
            KTMDA(n,i,:)=pc(:,i)*we(i);
            KTMDfA(n,i,:)=pc(:,i)*wef(i);
            KTMDrA(n,i,:)=pc(:,i)*wer(i);
        end
        n=n+1;
    end
    KTMDA=mean(squeeze(mean(KTMDA)));
    KTMDfA=mean(squeeze(mean(KTMDfA)));
    KTMDrA=mean(squeeze(mean(KTMDrA)));
    
    n=1;
    for s=[4 6 7]
        pc=squeeze(coe(8,s,:,:));
        we=mean(KTMDRPCA(1+20*(s-1):20+20*(s-1),:));
        wef=mean(KTMDRfPCA(1+20*(s-1):20+20*(s-1),:));
        wer=mean(KTMDRrPCA(1+20*(s-1):20+20*(s-1),:));
        for i=1:NPC
            KTMDRA(n,i,:)=pc(:,i)*we(i);
            KTMDRfA(n,i,:)=pc(:,i)*wef(i);
            KTMDRrA(n,i,:)=pc(:,i)*wer(i);
        end
        n=n+1;
    end
    KTMDRA=mean(squeeze(mean(KTMDRA)));
    KTMDRfA=mean(squeeze(mean(KTMDRfA)));
    KTMDRrA=mean(squeeze(mean(KTMDRrA))); 
    
    save results_FCchronosBlock_Monkey_Kin2.mat FowRev FowRevPCA NSUB coe lat ...
    WakePCA SleepPCA KTWPCA KTPCA KTRPCA KTMDWPCA KTMDPCA KTMDRPCA ...
    PFWPCA PFPCA PFRPCA MDWPCA MDPCA MDRPCA ...
    KTMDWA KTMDA KTMDRA ...
    KTMDWfA KTMDfA KTMDRfA ...
    KTMDWrA KTMDrA KTMDRrA ...
    Wakehier Sleephier KTWhier KThier KTRhier KTMDWhier KTMDhier KTMDRhier ...
    PFWhier PFhier PFRhier MDWhier MDhier MDRhier;
end
 
%%%%%%%

%% Rendering Su

if Su==1
    clear WakeA SleepA KTWA KTA KTRA KTMDWA KTMDA KTMDRA ...
    PFWA PFA PFRA MDWA MDA MDRA ...
    WakefA SleepfA KTWfA KTfA KTRfA KTMDWfA KTMDfA KTMDRfA ...
    PFWfA PFfA PFRfA MDWfA MDfA MDRfA ...
    WakerA SleeprA KTWrA KTrA KTRrA KTMDWrA KTMDrA KTMDRrA ...
    PFWrA PFrA PFRrA MDWrA MDrA MDRrA;   
    n=1;
    for s=[5 8 9]
        pc=squeeze(coe(6,s,:,:));
        we=mean(KTMDWPCA(1+20*(s-1):20+20*(s-1),:));
        wef=mean(KTMDWfPCA(1+20*(s-1):20+20*(s-1),:));
        wer=mean(KTMDWrPCA(1+20*(s-1):20+20*(s-1),:));
        for i=1:NPC
            KTMDWA(n,i,:)=pc(:,i)*we(i);
            KTMDWfA(n,i,:)=pc(:,i)*wef(i);
            KTMDWrA(n,i,:)=pc(:,i)*wer(i);
        end
        n=n+1;
    end
    KTMDWA=mean(squeeze(mean(KTMDWA)));
    KTMDWfA=mean(squeeze(mean(KTMDWfA)));
    KTMDWrA=mean(squeeze(mean(KTMDWrA)));
    
    n=1;
    for s=[5 8 9]
        pc=squeeze(coe(7,s,:,:));
        we=mean(KTMDPCA(1+20*(s-1):20+20*(s-1),:));
        wef=mean(KTMDfPCA(1+20*(s-1):20+20*(s-1),:));
        wer=mean(KTMDrPCA(1+20*(s-1):20+20*(s-1),:));
        for i=1:NPC
            KTMDA(n,i,:)=pc(:,i)*we(i);
            KTMDfA(n,i,:)=pc(:,i)*wef(i);
            KTMDrA(n,i,:)=pc(:,i)*wer(i);
        end
        n=n+1;
    end
    KTMDA=mean(squeeze(mean(KTMDA)));
    KTMDfA=mean(squeeze(mean(KTMDfA)));
    KTMDrA=mean(squeeze(mean(KTMDrA)));
    
    n=1;
    for s=[5 8 9]
        pc=squeeze(coe(8,s,:,:));
        we=mean(KTMDRPCA(1+20*(s-1):20+20*(s-1),:));
        wef=mean(KTMDRfPCA(1+20*(s-1):20+20*(s-1),:));
        wer=mean(KTMDRrPCA(1+20*(s-1):20+20*(s-1),:));
        for i=1:NPC
            KTMDRA(n,i,:)=pc(:,i)*we(i);
            KTMDRfA(n,i,:)=pc(:,i)*wef(i);
            KTMDRrA(n,i,:)=pc(:,i)*wer(i);
        end
        n=n+1;
    end
    KTMDRA=mean(squeeze(mean(KTMDRA)));
    KTMDRfA=mean(squeeze(mean(KTMDRfA)));
    KTMDRrA=mean(squeeze(mean(KTMDRrA))); 
    
    save results_FCchronosBlock_Monkey_Su.mat FowRev FowRevPCA NSUB coe lat ...
    WakePCA SleepPCA KTWPCA KTPCA KTRPCA KTMDWPCA KTMDPCA KTMDRPCA ...
    PFWPCA PFPCA PFRPCA MDWPCA MDPCA MDRPCA ...
    KTMDWA KTMDA KTMDRA ...
    KTMDWfA KTMDfA KTMDRfA ...
    KTMDWrA KTMDrA KTMDRrA ...
    Wakehier Sleephier KTWhier KThier KTRhier KTMDWhier KTMDhier KTMDRhier ...
    PFWhier PFhier PFRhier MDWhier MDhier MDRhier;
end
 
%%

%% Violin plots

THRLOW=5;
THRHIGH=90;

for sub=1:NSUB(1)
    Wake0(sub)=FowRev{1,sub};
end
Wake=rmoutliers(Wake0,'percentiles',[THRLOW THRHIGH]);
NSUB2(1)=length(Wake);
for sub=1:NSUB(2)
    Sleep0(sub)=FowRev{2,sub};
end
Sleep=rmoutliers(Sleep0,'percentiles',[THRLOW THRHIGH]);
NSUB2(2)=length(Sleep);
for sub=1:NSUB(3)
    KTW0(sub)=FowRev{3,sub};
end
KTW=rmoutliers(KTW0,'percentiles',[THRLOW THRHIGH]);
NSUB2(3)=length(KTW);
for sub=1:NSUB(4)
    KT0(sub)=FowRev{4,sub};
end
KT=rmoutliers(KT0,'percentiles',[THRLOW THRHIGH]);
NSUB2(4)=length(KT);
for sub=1:NSUB(5)
    KTR0(sub)=FowRev{5,sub};
end
KTR=rmoutliers(KTR0,'percentiles',[THRLOW THRHIGH]);
NSUB2(5)=length(KTR);
for sub=1:NSUB(6)
    KTMDW0(sub)=FowRev{6,sub};
end
KTMDW=rmoutliers(KTMDW0,'percentiles',[THRLOW THRHIGH]);
NSUB2(6)=length(KTMDW);
for sub=1:NSUB(7)
    KTMD0(sub)=FowRev{7,sub};
end
KTMD=rmoutliers(KTMD0,'percentiles',[THRLOW THRHIGH]);
NSUB2(7)=length(KTMD);
for sub=1:NSUB(8)
    KTMDR0(sub)=FowRev{8,sub};
end
KTMDR=rmoutliers(KTMDR0,'percentiles',[THRLOW THRHIGH]);
NSUB2(8)=length(KTMDR);
for sub=1:NSUB(9)
    PFW0(sub)=FowRev{9,sub};
end
PFW=rmoutliers(PFW0,'percentiles',[THRLOW THRHIGH]);
NSUB2(9)=length(PFW);
for sub=1:NSUB(10)
    PF0(sub)=FowRev{10,sub};
end
PF=rmoutliers(PF0,'percentiles',[THRLOW THRHIGH]);
NSUB2(10)=length(PF);
for sub=1:NSUB(11)
    PFR0(sub)=FowRev{11,sub};
end
PFR=rmoutliers(PFR0,'percentiles',[THRLOW THRHIGH]);
NSUB2(11)=length(PFR);
for sub=1:NSUB(12)
    MDW0(sub)=FowRev{12,sub};
end
MDW=rmoutliers(MDW0,'percentiles',[THRLOW THRHIGH]);
NSUB2(12)=length(MDW);
for sub=1:NSUB(13)
    MD0(sub)=FowRev{13,sub};
end
MD=rmoutliers(MD0,'percentiles',[THRLOW THRHIGH]);
NSUB2(13)=length(MD);
for sub=1:NSUB(14)
    MDR0(sub)=FowRev{14,sub};
end
MDR=rmoutliers(MDR0,'percentiles',[THRLOW THRHIGH]);
NSUB2(14)=length(MDR);



figure(1);
violinplot([Wake Sleep],[zeros(1,NSUB2(1)),ones(1,NSUB2(2))]);
title('Wake/Sleep');
savefig('1_violin_sleep.fig');
ranksum(Wake,Sleep)

figure(2);
violinplot([KTW' KT' KTR']);
title('Ketamine');
savefig('2_violin_KT.fig');
ranksum(KTW,KT)
ranksum(KT,KTR)

figure(3);
violinplot([KTMDW KTMD KTMDR],[zeros(1,NSUB2(6)),ones(1,NSUB2(7)),2*ones(1,NSUB2(8))]);
title('Ketamine-Medetodomine');
savefig('3_violin_KTMD.fig');
ranksum(KTMDW,KTMD)
ranksum(KTMD,KTMDR)

figure(4);
violinplot([PFW' PF' PFR']);
title('Propofol');
savefig('4_violin_PF.fig');
ranksum(PFW,PF)
ranksum(PF,PFR)

figure(5);
violinplot([MDW' MD' MDR']);
title('Medetodomine');
savefig('5_violin_MD.fig');
ranksum(MDW,MD)
ranksum(MD,MDR)


%% Time courses
WSvec=[Wake0(1:20) nan nan nan nan nan Sleep0(1:20)];
KTvec=[KTW0(1:20) nan nan nan nan nan KT0(1:20) nan nan nan nan nan KTR0(1:20)];
PFvec=[PFW0(21:40) nan nan nan nan nan PF0(21:40) nan nan nan nan nan PFR0(21:40)];
KTMDvec=[KTMDW0(61:80) nan nan nan nan nan KTMD0(61:80) nan nan nan nan nan KTMDR0(61:80)];
MDvec=[MDW0(1:20) nan nan nan nan nan MD0(1:20) nan nan nan nan nan MDR0(1:20)];

plot(WSvec);
savefig('1_ts_sleep.fig');
plot(KTvec);
savefig('2_ts_KT.fig');
plot(KTMDvec);
savefig('3_ts_KTMD.fig');
plot(PFvec);
savefig('4_ts_PF.fig');
plot(MDvec);
savefig('5_ts_MD.fig');

%% FC

task={'Wake';'Sleep';'KTW';'KT';'KTR';'KTMDW';'KTMD';'KTMDR';'PFW';'PF';'PFR';'MDW';'MD';'MDR'};

SES{1,:}=[1 2 3 8];  % Chibi Wake
SES{2,:}=[1 2 3 8];  % Chibi sleep
SES{3,:}=[1 4];    % Chibi KT
SES{4,:}=[1 4];
SES{5,:}=[1 4];
SES{6,:}=[1 2 3];  % George KTMD
SES{7,:}=[1 2 3];
SES{8,:}=[1 2 3];
SES{9,:}=[2 4];   % George PF
SES{10,:}=[2 4];
SES{11,:}=[2 4];
SES{12,:}=[1 3];   % Chibi MD
SES{13,:}=[1 3];
SES{14,:}=[1 3];



for xx=1:size(task,1)
    xx
    list=dir(['../../DataSet/MonkeyECoG/' task{xx} '/*.mat']);
    
    tss=[];
    for sub=SES{xx,:}
        sub
        load(list(sub).name);
        ts=TimeSeries;
        tss=[tss ts];
    end
    
    for seed=1:N
        tss(seed,:)=detrend(tss(seed,:)-nanmean(tss(seed,:)));
    end    
    X=tss';
    X = bsxfun(@minus,X,mean(X));
    [coe(xx,sub,:,:),pcats,lat(xx,sub,:)]=pca(X);
    tss2=pcats(:,1:NPC)';
    FCpca(xx,:,:)=cov(tss2');  %%corr(tss2(:,1:end-Tau)',tss2(:,1+Tau:end)');
    
    FC(xx,:,:)=corrcoef(tss');
end

%% FC
fc1=squeeze(FC(1,:,:));
fc2=squeeze(FC(2,:,:));
subplot(1,3,1);
imagesc(fc1);
axis square;
subplot(1,3,2);
imagesc(fc2);
axis square;
subplot(1,3,3);
scatter(fc1(:),fc2(:),'.');
axis square;
savefig('1_FCwakesleep.fig');

fc1=squeeze(FC(3,:,:));
fc2=squeeze(FC(4,:,:));
fc3=squeeze(FC(5,:,:));
subplot(1,4,1);
imagesc(fc1);
axis square;
subplot(1,4,2);
imagesc(fc2);
axis square;
subplot(1,4,3);
imagesc(fc3);
axis square;
subplot(1,4,4);
scatter(fc1(:),fc2(:),'.');
axis square;
savefig('2_FC_KT.fig');

fc1=squeeze(FC(6,:,:));
fc2=squeeze(FC(7,:,:));
fc3=squeeze(FC(8,:,:));
subplot(1,4,1);
imagesc(fc1);
axis square;
subplot(1,4,2);
imagesc(fc2);
axis square;
subplot(1,4,3);
imagesc(fc3);
axis square;
subplot(1,4,4);
scatter(fc1(:),fc2(:),'.');
axis square;
savefig('3_FC_KTMD.fig');

fc1=squeeze(FC(9,:,:));
fc2=squeeze(FC(10,:,:));
fc3=squeeze(FC(11,:,:));
subplot(1,4,1);
imagesc(fc1);
axis square;
subplot(1,4,2);
imagesc(fc2);
axis square;
subplot(1,4,3);
imagesc(fc3);
axis square;
subplot(1,4,4);
scatter(fc1(:),fc2(:),'.');
axis square;
savefig('4_FC_PF.fig');

fc1=squeeze(FC(12,:,:));
fc2=squeeze(FC(13,:,:));
fc3=squeeze(FC(14,:,:));
subplot(1,4,1);
imagesc(fc1);
axis square;
subplot(1,4,2);
imagesc(fc2);
axis square;
subplot(1,4,3);
imagesc(fc3);
axis square;
subplot(1,4,4);
scatter(fc1(:),fc2(:),'.');
axis square;
savefig('5_FC_MD.fig');

%% FC PCA
figure;
subplot(1,5,1)
fc1=squeeze(FCpca(1,:,:));
fc2=squeeze(FCpca(2,:,:));
fc1=diag(fc1);
fc2=diag(fc2);
scatter(fc1,fc2,'.');
maxf=max(fc1(1),fc2(1));
xlim([0 maxf]);
ylim([0 maxf]);
axis square;

subplot(1,5,5)
fc1=squeeze(FCpca(3,:,:));
fc2=squeeze(FCpca(4,:,:));
fc1=diag(fc1);
fc2=diag(fc2);
scatter(fc1,fc2,'.');
maxf=max(fc1(1),fc2(1));
xlim([0 maxf]);
ylim([0 maxf]);
axis square;

subplot(1,5,4)
fc1=squeeze(FCpca(6,:,:));
fc2=squeeze(FCpca(7,:,:));
fc1=diag(fc1);
fc2=diag(fc2);
scatter(fc1,fc2,'.');
maxf=max(fc1(1),fc2(1));
xlim([0 maxf]);
ylim([0 maxf]);
axis square;

subplot(1,5,2)
fc1=squeeze(FCpca(9,:,:));
fc2=squeeze(FCpca(10,:,:));
fc1=diag(fc1);
fc2=diag(fc2);
scatter(fc1,fc2,'.');
maxf=max(fc1(1),fc2(1));
xlim([0 maxf]);
ylim([0 maxf]);
axis square;

subplot(1,5,3)
fc1=squeeze(FCpca(12,:,:));
fc2=squeeze(FCpca(13,:,:));
fc1=diag(fc1);
fc2=diag(fc2);
scatter(fc1,fc2,'.');
maxf=max(fc1(1),fc2(1));
xlim([0 maxf]);
ylim([0 maxf]);
axis square;
savefig('FigSupl_FCpca.fig');

%% Violin plots Hierarchy

figure(1)
Wakehier2=rmoutliers(Wakehier,'percentiles',[THRLOW THRHIGH]);
Sleephier2=rmoutliers(Sleephier,'percentiles',[THRLOW THRHIGH]);
violinplot([Wakehier2 Sleephier2],[zeros(1,length(Wakehier2)) ones(1,length(Sleephier2))]);
ranksum(Wakehier,Sleephier)
title('Wake/Sleep');
savefig('1_violin_sleep_hier.fig');

figure(2)
KTWhier2=rmoutliers(KTWhier,'percentiles',[THRLOW THRHIGH]);
KThier2=rmoutliers(KThier,'percentiles',[THRLOW THRHIGH]);
KTRhier2=rmoutliers(KTRhier,'percentiles',[THRLOW THRHIGH]);
violinplot([KTWhier2 KThier2 KTRhier2],[zeros(1,length(KTWhier2)) ones(1,length(KThier2)) 2*ones(1,length(KTRhier2))]);
ranksum(KTWhier,KThier)
ranksum(KThier,KTRhier)
title('Ketamine');
savefig('2_violin_KT_hier.fig');

figure(3)
KTMDWhier2=rmoutliers(KTMDWhier,'percentiles',[THRLOW THRHIGH]);
KTMDhier2=rmoutliers(KTMDhier,'percentiles',[THRLOW THRHIGH]);
KTMDRhier2=rmoutliers(KTMDRhier,'percentiles',[THRLOW THRHIGH]);
violinplot([KTMDWhier2 KTMDhier2 KTMDRhier2],[zeros(1,length(KTMDWhier2)) ones(1,length(KTMDhier2)) 2*ones(1,length(KTMDRhier2))]);
ranksum(KTMDWhier,KTMDhier)
ranksum(KTMDhier,KTMDRhier)
title('KTMD');
savefig('3_violin_KTMD_hier.fig');

figure(4)
PFWhier2=rmoutliers(PFWhier,'percentiles',[THRLOW THRHIGH]);
PFhier2=rmoutliers(PFhier,'percentiles',[THRLOW THRHIGH]);
PFRhier2=rmoutliers(PFRhier,'percentiles',[THRLOW THRHIGH]);
violinplot([PFWhier2 PFhier2 PFRhier2],[zeros(1,length(PFWhier2)) ones(1,length(PFhier2)) 2*ones(1,length(PFRhier2))]);
ranksum(PFWhier,PFhier)
ranksum(PFhier,PFRhier)
title('PF');
savefig('4_violin_PF_hier.fig');

figure(5)
MDWhier2=rmoutliers(MDWhier,'percentiles',[THRLOW THRHIGH]);
MDhier2=rmoutliers(MDhier,'percentiles',[THRLOW THRHIGH]);
MDRhier2=rmoutliers(MDRhier,'percentiles',[THRLOW THRHIGH]);
violinplot([MDWhier2 MDhier2 MDRhier2],[zeros(1,length(MDWhier2)) ones(1,length(MDhier2)) 2*ones(1,length(MDRhier2))]);
ranksum(MDWhier,MDhier)
ranksum(MDhier,MDRhier)
title('MD');
savefig('5_violin_MD_hier.fig');
