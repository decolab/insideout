
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
    nn=1;
    for sub=SES{xx,:}
        load(list(sub).name);
        ts=TimeSeries;
        FC(xx,nn)=mean(mean(corrcoef(ts')));
        nn=nn+1;
    end
end

for i=1:14
    for j=i+1:14
        p=ranksum(FC(i,:),FC(j,:));
        if p<0.05
            i
            j
        end
    end
end
