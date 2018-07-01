%%Use this script to transform the matlab data into the format used in the R and JAGS script
%% Read in data from files
clear
root='....'; %where your data is stored
% select experiment
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
        subs   = dir([folder.name '\Data\Subject*']);                    %subjects' data
end

if length(subs)==0, error('no data found'); end
for i=1:length(subs);   % for each subject folder
    files = dir([folder.name '\Data\' subs(i).name '\']);            % find the file
    fname = [folder.name '\Data\' subs(i).name '\' subs(i).name '.mat']; % bild filename
    tmp = load( fname,'result' );                               % load file
    fprintf('loaded %s: %g\n', fname, length(tmp.result.data));
    if i>1                % not all subject data files are the same format...
        [tmp.result.data alld ] = ensureStructsAssignable(tmp.result.data, alld);
    end                   % so ensure they match
    alld(i,:) = tmp.result.data;                                % store data
end
fprintf('loaded %g subjects, %g trials\n', size(alld));
%%
arrow = catStruct(1,alld,'showProbeInitially'); % was it an arrow probe?
switch EXPT
    
    case 1
        cond  = catStruct(1,alld,'nTargets');   % condition: 3 or 5 targets
        cond  = arrow + (cond-1)*2 + 1;         % probe type AND num targets
        cond2 = (catStruct(1,alld,'nTargets')-1)/2 + 4;
        NC    = 2; % two conditions, after the probe-type effect
        maxitems = 5;
        conditionNames = {'3 circ', '3 arr', '5 circ', '5 arr'};
    case 2
        cond  = catStruct(1,alld,'drawMask'); % condition: mask present?
        cond  = arrow + cond*2 + 1;           % probe type AND mask present
        cond2 = catStruct(1,alld,'drawMask') + 5;
        NC    = 2;
        maxitems = 3;
        conditionNames = {'nomask circ','nomask arr','mask circ','mask arr'};
    case 3
        seq   = permute(catStruct(3,alld,'sequence'),[4,2,3,1]);
        ti    = sq(sum( (seq==1) .* repmat(1:4,[288,1,14]) , 2 ))'; % target index
        cond  = arrow + (ti-1)*2 + 1;            % probe type AND seq pos
        cond2 = ti+8;
        NC    = 4; % conditions, other than probe type, = 4 sequential positions
        maxitems = 4;
        conditionNames = {'1 circ','1 arr', '2 circ','2 arr', '3 circ', '3 arr', '4 circ', '4 arr'};
        
        
end

%% Set up parameters for Bayes Hierarchical Model
memset=NaN(size(alld,1)*size(alld,2),maxitems);
k=1;
for subject=1:size(alld,1) %for all subjects
    for trial=1:size(alld,2) %for all trials
        memset(k,1:length(alld(subject,trial).angles))=alld(subject,trial).angles; %presented stimuli angles (first column is target stimulus)
        x(k,:)=alld(subject,trial).response;         %response angle
        id(k,:)=subject;                             %participant id
        condition(k,:)=cond(subject,trial);          %condition number
        setsize(k,:) = alld(subject,trial).nTargets; %setsize
        k=k+1;
    end
end

if EXPT==1
    condition(condition==5)  = 1; %rename conditions to 1-4
    condition(condition==6)  = 2; 
    condition(condition==9)  = 3; 
    condition(condition==10) = 4; 
end

%save all of those variables as .dat files

save([root '/Experiment ' num2str(EXPT) '_m.dat'],'memset','-ascii')
save([root '/Experiment ' num2str(EXPT) '_x.dat'],'x','-ascii')
save([root '/Experiment ' num2str(EXPT) '_id.dat'],'id','-ascii')
save([root '/Experiment ' num2str(EXPT) '_condition.dat'],'condition','-ascii')
save([root '/Experiment ' num2str(EXPT) '_setsize.dat'],'setsize','-ascii')


