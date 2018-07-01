%Recall cues interfere with  retrieval from visuospatial working memory
%% Read in data from files
clear
root='C:\...'; %where you stored the data
% select experimen
EXPT = 1 %which experiment do you want to analyse?
switch EXPT  % which folder and files?
    case 1     % 3/5 items simultaneous, probe by colour
        folder.name = [root 'Experiment 1: Simultaneous Presentation']; %folder with files
        subs   = dir([folder.name '\Data\Subject*']);                   %subjects' data
    case 2     % masking/no masking, probe by colour
        folder.name = [root 'Experiment 2: Masking'];                   %folder with files
        subs   = dir([folder.name '\Data\Subject*']);                   %subjects' data
    case 3     % sequential presentation, probe by colour
        folder.name = [root 'Experiment 3: Sequential Presentation'];   %folder with files
        subs   = dir([folder.name '\Data\Subject*'])                    %subjects' data
end

if length(subs)==0, error('no data found'); end                         %if subject data is empty
for i=1:length(subs);                                                   % for each subject folder
    files = dir([folder.name '\Data\' subs(i).name '\']);               % find the file
    fname = [folder.name '\Data\' subs(i).name '\' subs(i).name '.mat'];% bild filename
    tmp = load( fname,'result' );                                       % load file
    fprintf('loaded %s: %g\n', fname, length(tmp.result.data));
    if i>1                % not all subject data files are the same format...
        [tmp.result.data alld ] = ensureStructsAssignable(tmp.result.data, alld);
    end                   % so ensure they match
    alld(i,:) = tmp.result.data;                                        % store data
end
fprintf('loaded %g subjects, %g trials\n', size(alld));

%% calculate response, target, probe, RT
arrow = catStruct(1,alld,'showProbeInitially'); % was it an arrow probe?
switch EXPT
    case 1
        cond  = catStruct(1,alld,'nTargets');   % condition: 3 or 5 targets
        cond  = arrow + (cond-1)*2 + 1;         % probe type AND num targets
        NC    = 2; % two conditions, after the probe-type effect
        conditionNames = {'3 circ', '3 arr', '5 circ', '5 arr'};
    case 2
        cond = catStruct(1,alld,'drawMask'); % condition: mask present?
        cond = arrow + cond*2 + 1;           % probe type AND mask present
        NC   = 2;
        conditionNames = {'nomask circ','nomask arr','mask circ','mask arr'};
    case 3
        seq  = permute(catStruct(3,alld,'sequence'),[4,2,3,1]);
        ti   = sq(sum( (seq==1) .* repmat(1:4,[288,1,15]) , 2 ))'; % target index
        cond = arrow + (ti-1)*2 + 1;            % probe type AND seq pos
        NC   = 4; % conditions, other than probe type, = 4 sequential positions
        conditionNames = {'1 circ','1 arr', '2 circ','2 arr', '3 circ', '3 arr', '4 circ', '4 arr'};
end

%get initiation starttime via trajectory which stored mouse movement over
%time. It's stored in the row and coloumn of the field trajectory 
for i=1:size(alld,1) 
    for j=1:size(alld,2)
        if ~isempty(alld(i,j).trajectory)
            starttime(i,j)=alld(i,j).trajectory(1,1);
        else
            starttime(i,j)=NaN;
        end
    end
end

r         = catStruct(1,alld,'response');            % r = response angle
t         = catStruct(1,alld,'targetAngle');         % t = target angle
p         = catStruct(1,alld,'responseOffset');      % p = target-to-initial-probe distance
rt        = catStruct(1, alld,'endChoice') - catStruct(1, alld,'startProbe'); % reaction time
startresp = starttime - catStruct(1, alld,'startProbe'); %time to initiate response

%% Plot Histograms of Reaction Time

hist(rt(:),500)
xlim([0 8])
vline(mean2(rt)+3*std2(rt),'r','Mean + 3*STD') %vertical line at x = reaction time mean plus 3*standard deviation

%% Bays Correction

%creating T and X

for subnum=1:length(subs) %for all subjects
    n=1;
    for condition=unique(cond)' %for all conditions
        X{subnum}{n}=r(subnum,cond(subnum,:)==condition & rt(subnum,:)<5)';    % create X, a (nx1)-vector of responses
        b=find(cond(subnum,:)==condition & rt(subnum,:)<5);                    % restricted to trials with reaction times under 5s
        for i=1:length(b)
            T{subnum}{n}(:,i)=(alld(subnum,b(i)).angles);  % create T, an (nxm) matrix of stimulus values ...
        end                                                % (the first line contains target values, ...
        n=n+1;                                             % the remaining lines non-target values))
        excl(subnum,condition) = sum(rt(subnum,cond(subnum,:)==condition)>5); %store which trials have been excluded
    end
end
excl=excl(:,unique(cond)');
clear n i condition b subnum


%% Mixture model analysis of analogue report data
%correcting by multiplying target frequence and absolute error

ae    = abs(wrap(r-t));   %absolute error
bad   = rt>5;             % exclude overly long trials
ae    = ae+bool2nan( bad ); % eg attentional lapses

condnum=unique(cond)';
for subnum=1:length(subs)                  %for all subjects
    for condition=1:length(conditionNames) %for all conditions
        [B LL W]=CO16_fit((X{subnum}{condition}-pi),(T{subnum}{condition}(1,:)'-pi),(T{subnum}{condition}(2:end,:)'-pi)); % CO16 function by Paul Bays (source: http://www.paulbays.com/code.php)
        [P(subnum,condition) B_error(subnum,condition)]=JV10_error((X{subnum}{condition}-pi),(T{subnum}{condition}(1,:)'-pi)); % JV10_error function by Paul Bays (source: http://www.paulbays.com/code.php)
        
        statistics(subnum,condition)=nanmean(ae(subnum,find(cond(subnum,:)==condnum(condition) & rt(subnum,:)<5)));             %mean absolute error
        startresp_stat(subnum,condition)=nanmean(startresp(subnum,find(cond(subnum,:)==condnum(condition) & rt(subnum,:)<5)));  %mean initiation time
        rt_stat(subnum,condition)=mean(rt(subnum,find(cond(subnum,:)==condnum(condition) & rt(subnum,:)<5)));                   %mean reaction time
        
        pT_stat(subnum,condition)=B(2);  %mean target response probability
        pNT_stat(subnum,condition)=B(3); %mean non-target response probability
        pU_stat(subnum,condition)=B(4);  %mean uniform response probability
        K_stat(subnum,condition)=B(1);   %mean distribution concentration parameter
        clear B LL %A holds the
    end
end
B=cat(3, pT_stat, pNT_stat, pU_stat, K_stat);
% B( subject , condition, T/NT/U/K )

clear condition condnum subnum
%% scatter plots of effects for experiment 1

variables = {'P_{target}','P_{nontarget}','P_{uniform}','\kappa'};
for i=1:length(variables)
    subplot(3,4,i) %horizontal: probabilities and concentration parameter
                   %vertical: 1. row: 3 arrow vs dot
                             %2. row: 5 arrow vs dot
                             %3. row: 5 items (dot-arrow) vs 3 items (dot-arrow)
    
    scatter(B(:,1,i),B(:,2,i),'filled'); title(variables{i});

    xlabel 'three items, dot'; ylabel 'three items, arrow'
    
    %adjust axis limits
    
    if i==1
        xlim([0 1]);
        ylim([0 1]);
    elseif i==4
        xlim([0 30]);
        ylim([0 30]);
    else
        xlim([0 0.4]);
        ylim([0 0.4]);
    end
    
    %plot diagonal
    hold on; plot(xlim,xlim,':','LineWidth',0.7); hold off
    %plot errorbars around the mean
    hold on; errorbar(mean(B(:,1,i)),mean(B(:,2,i)),...
        std(B(:,2,i))/2,std(B(:,2,i))/2,...
        std(B(:,1,i))/2,std(B(:,1,i))/2,'o','LineWidth',1); hold off
    
    subplot(3,4,i+4)
    scatter(B(:,3,i),B(:,4,i),'filled'); title(variables{i});
    
    xlabel 'five items, dot'; ylabel 'five items, arrow'
    
    %adjust axis limits    
    
    if i<4
        xlim([0 1]);
        ylim([0 1]);
    else
        xlim([0 30]);
        ylim([0 30]);
    end

    %plot diagonal
    hold on; plot(xlim,xlim,':','LineWidth',0.7); hold off
    %plot errorbars around mean
    hold on; errorbar(mean(B(:,3,i)),mean(B(:,4,i)),...
        std(B(:,4,i))/2,std(B(:,4,i))/2,...
        std(B(:,3,i))/2,std(B(:,3,i))/2,'o','LineWidth',1); hold off
    
    subplot(3,4,i+8)
    scatter(B(:,1,i)-B(:,2,i),B(:,3,i)-B(:,4,i),'filled');
    
    %adjust axis limits
    if i<4
        xlim([-0.3 0.3]);
        ylim([-0.6 0.6]);
    else
        xlim([-10 10]);
        ylim([-10 10]);
    end
    
    %plot diagonal
    hold on; plot(xlim,[0 0],':','LineWidth',0.7); plot([0 0],ylim,':','LineWidth',0.7);hold off
    %plot errorbars around mean
    hold on; errorbar(mean(B(:,1,i)-B(:,2,i)),mean(B(:,3,i)-B(:,4,i)),...
        std(B(:,3,i)-B(:,4,i))/2,std(B(:,3,i)-B(:,4,i))/2,...
        std(B(:,1,i)-B(:,2,i))/2,std(B(:,1,i)-B(:,2,i))/2,'o','Color','red','LineWidth',1); hold off
    title(variables{i});
    
    xlabel 'three items (dot-arr)'; ylabel 'five items (dot-arr)'
    
end

%% Plot

ae    = abs(wrap(r-t));     % absolute error
NS    = size(ae,1);         % num subjects
% exclude trials?
bad   = rt>5;             % exclude overly long trials
ae    = ae+bool2nan( bad ); % eg attentional lapses
% group abs error by condition, using groupMeans(ae, cond)
ae   = nancat([2,3],nancat([1,2],groupMeans(ae,2, cond,@(x){x}))); ae=ae{1};
% now we have abs err (TRIAL, SUBJECT, COND)
if 0 % plot all conditions in a row (flat)?
    errorBarPlot(sq(nanmean(ae)));
elseif 1 % plot conditions as dimensions -- As Lines?
    errorBarPlot(permute(reshape(sq(nanmean(ae)),NS,2,NC),[1,3,2]));
else  % conditions as dimensions -- As bars?
    errorBarPlot(permute(reshape(sq(nanmean(ae)),NS,2,NC),[1,3,2]),'type','bar');
end
set(gca,'xtick',1:length(conditionNames),'xticklabel',conditionNames);;
legend('Circle','Arrow')
title 'Mean abs error'; drawnow

%% Plot Bays: Non-Target Response
figure
errorBarPlot(permute(reshape(pNT_stat,NS,2,NC),[1,3,2]));
set(gca,'xtick',1:length(conditionNames),'xticklabel',conditionNames);;
legend('Circle','Arrow')
title 'Non-Target Responses'; drawnow

%% Plot Bays: Uniform Response
figure
errorBarPlot(permute(reshape(pU_stat,NS,2,NC),[1,3,2]));
set(gca,'xtick',1:length(conditionNames),'xticklabel',conditionNames);;
legend('Circle','Arrow')
title 'Uniform Responses'; drawnow

%% Plot Bays: Concentration Parameter
figure
errorBarPlot(permute(reshape(K_stat,NS,2,NC),[1,3,2]));
set(gca,'xtick',1:length(conditionNames),'xticklabel',conditionNames);;
legend('Circle','Arrow')
title 'concentration parameter'; drawnow
%% Plot RT
NS    = size(startresp_stat,1);
figure
errorBarPlot(permute(reshape(startresp_stat,NS,2,NC),[1,3,2]));
legend('Circle','Arrow')
title 'initiation time'; drawnow

figure
errorBarPlot(permute(reshape(rt_stat,NS,2,NC),[1,3,2]));
legend('Circle','Arrow')
title 'finishing time'; drawnow

%% plot quantile contours of error distribution, as a function of probe angle
colourMap([jet(40);flipud(jet(40))]);
sel=bool2nan(~arrow);clf
[h,xbin,ybin]=plotBinsQuantiled( wrap(p+sel)', wrap(r-t+sel)',160,'Bins',2) ;
hold on; plot(xlim,[0 0],':');
%% Empirical (Kaplan-Meier) cumulative distribution function for experiment 3

err    = (wrap(r-t));     % absolute error
NS    = size(err,1);         % num subjects
% exclude trials?
bad   = rt>5;             % exclude overly long trials
err    = err+bool2nan( bad ); % eg attentional lapses
% group abs error by condition, using groupMeans(ae, cond)
err   = nancat([2,3],nancat([1,2],groupMeans(err,2, cond,@(x){x}))); err=err{1};
%uniform distribution in the circular space
uniformfun=@(x) (x+pi)/(2*pi);

%use ks-tests to test whether the items in different sequential positions
%are unequal to the uniform distribution
for i=1:size(err,3)
    singlerr{i} = err(:,:,i);
    [f{i},x{i}]=ecdf(singlerr{i}(:));
    [h_ks{i},p_ks{i}]=kstest(singlerr{i}(:),'CDF',[singlerr{i}(:),cdf('unif',singlerr{i}(:),-pi,+pi)]);
end

hold on
plot(x{1},f{1},'--','Linewidth',2) %plot circle condition of first sequentially presented item
plot(x{2},f{2},'r','Linewidth',2)  %plot arrow condition of first sequentially presented item
fplot(uniformfun, [-pi pi], 'k','Linewidth',2); %plot uniform distribution
hold off
title 'ecdf of the first sequential position'; drawnow
legend('circle','arrow','uniform')


%% Context effect for experiment 1
NS    = size(arrow,1);
context=[NaN(15,1) arrow]+2*[arrow NaN(15,1)]; %add a matrix shifted by 1 to the original matrix
context(:,end)=[];
%(+3) - arrow, arrow
%(+2) - circle, arrow
%(+1) - arrow, circle
%(+0) - circle, circle
nTargets         = catStruct(1,alld,'nTargets');  % extract number of targets

for subjects=1:15
    contmat(subjects,1)=nanmean(ae(subjects,context(subjects,:)==3 & nTargets(subjects,:)==3)'); %three items, arrow preceded by arrow
    contmat(subjects,2)=nanmean(ae(subjects,context(subjects,:)==2 & nTargets(subjects,:)==3)'); %three items, arrow preceded by circle
    contmat(subjects,3)=nanmean(ae(subjects,context(subjects,:)==1 & nTargets(subjects,:)==3)'); %three items, circle preceded by arrow
    contmat(subjects,4)=nanmean(ae(subjects,context(subjects,:)==0 & nTargets(subjects,:)==3)'); %three items, circle preceded by circle
    contmat(subjects,5)=nanmean(ae(subjects,context(subjects,:)==3 & nTargets(subjects,:)==5)'); %five items, arrow preceded by arrow
    contmat(subjects,6)=nanmean(ae(subjects,context(subjects,:)==2 & nTargets(subjects,:)==5)'); %five items, arrow preceded by circle
    contmat(subjects,7)=nanmean(ae(subjects,context(subjects,:)==1 & nTargets(subjects,:)==5)'); %five items, circle preceded by arrow
    contmat(subjects,8)=nanmean(ae(subjects,context(subjects,:)==0 & nTargets(subjects,:)==5)'); %five items, circle preceded by circle
end

errorBarPlot(permute(reshape(contmat,15,2,4),[1,3,2])); %plot context graph
xticklabels({'3, arrow','3, circle','5, arrow','5, circle'})
legend('circle','arrow','uniform')
legend('Preceded by Arrow','Preceded by Circle')


