clear all;
path2=[ '../TENET/'];
addpath(path2);

tasks={'REST1';'EMOTION';'GAMBLING';'WM';'LANGUAGE';'MOTOR';'RELATIONAL';'SOCIAL'};

N=80;
NSUB=990;

NLAG=6;

FowRev=zeros(size(tasks,1),NLAG,NSUB);
AsymFow=zeros(size(tasks,1),NLAG,NSUB);
AsymRev=zeros(size(tasks,1),NLAG,NSUB);

for xx=1:size(tasks,1)
    xx
    load(['hcp1003_' tasks{xx} '_LR_dbs80.mat']);
    
    nsub=1;
    for sub=1:1003
        if isstruct(subject{sub})
            if size(subject{sub}.dbs80ts,2)>175 % one subject has only in EMOTION
                subject{nsub}=subject{sub};
                nsub=nsub+1;
            end
        end
        
    end
%     for sub=1:nsub-1
%         subject{sub}.dbs80ts=subject{sub}.dbs80ts(:,1:176);
%     end
    
    
    for sub=1:NSUB
        sub
        ts=subject{sub}.dbs80ts;
        Tm=size(subject{sub}.dbs80ts,2);
        for Tau=1:NLAG
            [FCtau_foward(Tau,:,:) pctauf]=corr(ts(:,1:Tm-Tau)',ts(:,1+Tau:Tm)');
            [FCtau_reversal(Tau,:,:) pctaur]=corr(ts(:,Tm:-1:Tau+1)',ts(:,Tm-Tau:-1:1)');
            FCtf=squeeze(FCtau_foward(Tau,:,:));
            FCtr=squeeze(FCtau_reversal(Tau,:,:));
            %             pctauf(pctauf>0.01)=2;
            %             pctauf(pctauf<=0.01)=1;
            %             pctauf(pctauf==2)=0;
            %             pctaur(pctaur>0.01)=2;
            %             pctaur(pctaur<=0.01)=1;
            %             pctaur(pctaur==2)=0;
            %             Itauf=-0.5*log(1-FCtf.*FCtf);
            %             Itaur=-0.5*log(1-FCtr.*FCtr);
            %             mutinf=pctauf(:).*pctaur(:).*(Itauf(:)-Itaur(:)).^2;
            Itauf=-0.5*log(1-FCtf.*FCtf);
            Itaur=-0.5*log(1-FCtr.*FCtr);
            Reference=((Itauf(:)-Itaur(:)).^2)';
            index=find(Reference>quantile(Reference,0.95));
            FowRev(xx,Tau,sub)=nanmean(Reference(index));
            AsymFow(xx,Tau,sub)=mean(mean(abs(Itauf-Itauf')));
            AsymRev(xx,Tau,sub)=mean(mean(abs(Itaur-Itaur')));
        end
    end
end

for xx=2:size(tasks,1)
    [aux index2(xx-1)]=max(mean(squeeze(FowRev(xx,:,:)),2));
end
Tauwinner=round(mean(index2))

for xx=2:size(tasks,1)
    ptau(xx-1)=ranksum(squeeze(FowRev(1,Tauwinner,:)),squeeze(FowRev(xx,Tauwinner,:)));
end

figure(1);
T=Tauwinner;
violinplot([squeeze(FowRev(1,T,:)) squeeze(FowRev(4,T,:)) squeeze(FowRev(5,T,:)) squeeze(FowRev(6,T,:)) squeeze(FowRev(3,T,:)) squeeze(FowRev(2,T,:)) squeeze(FowRev(7,T,:)) squeeze(FowRev(8,T,:))]);
ranksum(squeeze(FowRev(1,T,:)),squeeze(FowRev(4,T,:)))

save results_FCchronos_REST1_Block.mat FowRev AsymFow AsymRev ptau Tauwinner;

