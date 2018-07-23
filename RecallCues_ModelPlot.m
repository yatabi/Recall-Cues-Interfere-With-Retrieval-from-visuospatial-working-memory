clear all
close all

Labels={'Three Items','Five Items';'No Mask','Mask'};
maxprecision=[22 22 22];
NC=[2 2 4];
root='...';
EXPT=3;
load([root 'EXPT_' num2str(EXPT) '_effects_1.mat']);

if EXPT~=3

    ptarget(:,2)=samples(:,3);
    ptarget(:,1)=(exp(samples(:,38)).*samples(:,3)./((exp(samples(:,38))-1).*samples(:,3)+1));
    ptarget(:,4)=(exp(samples(:,106)).*samples(:,3)./((exp(samples(:,106))-1).*samples(:,3)+1));
    ptarget(:,3)=(exp(samples(:,72)).*samples(:,3)./((exp(samples(:,72))-1).*samples(:,3)+1));
    
    precision(:,2)=samples(:,1);
    precision(:,1)=samples(:,1)+samples(:,36);
    precision(:,4)=samples(:,1)+samples(:,104);
    precision(:,3)=samples(:,1)+samples(:,36)+samples(:,104)+samples(:,70);
    
else
    ptarget(:,2)=samples(:,3);
    ptarget(:,1)=(exp(samples(:,36)).*samples(:,3)./((exp(samples(:,36))-1).*samples(:,3)+1));
    ptarget(:,4)=(exp(samples(:,4539)).*samples(:,3)./((exp(samples(:,4539))-1).*samples(:,3)+1));
    ptarget(:,3)=(exp(samples(:,72)).*samples(:,3)./((exp(samples(:,72))-1).*samples(:,3)+1));
    ptarget(:,6)=(exp(samples(:,4573)).*samples(:,3)./((exp(samples(:,4573))-1).*samples(:,3)+1));
    ptarget(:,5)=(exp(samples(:,106)).*samples(:,3)./((exp(samples(:,106))-1).*samples(:,3)+1));
    ptarget(:,8)=(exp(samples(:,4607)).*samples(:,3)./((exp(samples(:,4607))-1).*samples(:,3)+1));
    ptarget(:,7)=(exp(samples(:,140)).*samples(:,3)./((exp(samples(:,140))-1).*samples(:,3)+1));
    
    precision(:,2)=samples(:,1);
    precision(:,1)=samples(:,1)+samples(:,36);
    precision(:,4)=samples(:,1)+samples(:,4537);
    precision(:,3)=samples(:,1)+samples(:,36)+samples(:,4537)+samples(:,70);
    precision(:,6)=samples(:,1)+samples(:,4571);
    precision(:,5)=samples(:,1)+samples(:,36)+samples(:,4571)+samples(:,104);
    precision(:,8)=samples(:,1)+samples(:,4605);
    precision(:,7)=samples(:,1)+samples(:,36)+samples(:,4605)+samples(:,138);
end


nontarget=(1-ptarget).*samples(:,2);
ptarget=samples(:,2).*ptarget;
uniform=repmat(1-samples(:,2),1,2*NC);
Legends={'Circle','Arrow'};

for i=1:2*NC(EXPT)
    HDIlim(i,:)=HDIofMCMC(ptarget(:,i),0.95);
end

subplot(2,2,1);
hold on
if EXPT~=3
errorbar([1 2],[mean(ptarget(:,1)) mean(ptarget(:,3))],[mean(ptarget(:,1))-HDIlim(1,1) mean(ptarget(:,3))-HDIlim(3,1)],[HDIlim(1,2)-mean(ptarget(:,1)) HDIlim(3,2)-mean(ptarget(:,3))])
errorbar([1 2],[mean(ptarget(:,2)) mean(ptarget(:,4))],[mean(ptarget(:,2))-HDIlim(2,1) mean(ptarget(:,4))-HDIlim(4,1)],[HDIlim(2,2)-mean(ptarget(:,2)) HDIlim(4,2)-mean(ptarget(:,4))])
else
errorbar([1 2 3 4],[mean(ptarget(:,1)) mean(ptarget(:,3)) mean(ptarget(:,5)) mean(ptarget(:,7))],[mean(ptarget(:,1))-HDIlim(1,1) mean(ptarget(:,3))-HDIlim(3,1) mean(ptarget(:,5))-HDIlim(5,1) mean(ptarget(:,7))-HDIlim(7,1)],[HDIlim(1,2)-mean(ptarget(:,1)) HDIlim(3,2)-mean(ptarget(:,3)) HDIlim(5,2)-mean(ptarget(:,5)) HDIlim(7,2)-mean(ptarget(:,7))])
errorbar([1 2 3 4],[mean(ptarget(:,2)) mean(ptarget(:,4)) mean(ptarget(:,6)) mean(ptarget(:,8))],[mean(ptarget(:,2))-HDIlim(2,1) mean(ptarget(:,4))-HDIlim(4,1) mean(ptarget(:,6))-HDIlim(6,1) mean(ptarget(:,8))-HDIlim(8,1)],[HDIlim(2,2)-mean(ptarget(:,2)) HDIlim(4,2)-mean(ptarget(:,4)) HDIlim(6,2)-mean(ptarget(:,6)) HDIlim(8,2)-mean(ptarget(:,8))])
end
hold off
ylim([0 1]);
xlim([0.8 (NC(EXPT)+0.2)]);
xticks([1:NC(EXPT)]);
legend(Legends);
if EXPT~=3
    xticklabels(Labels(EXPT,:));
end
title('Target Response');

for i=1:2*NC(EXPT)
  HDIlim(i,:)=HDIofMCMC(nontarget(:,i),0.95);
end

subplot(2,2,2);
hold on
if EXPT~=3
errorbar([1 2],[mean(nontarget(:,1)) mean(nontarget(:,3))],[mean(nontarget(:,1))-HDIlim(1,1) mean(nontarget(:,3))-HDIlim(3,1)],[HDIlim(1,2)-mean(nontarget(:,1)) HDIlim(3,2)-mean(nontarget(:,3))])
errorbar([1 2],[mean(nontarget(:,2)) mean(nontarget(:,4))],[mean(nontarget(:,2))-HDIlim(2,1) mean(nontarget(:,4))-HDIlim(4,1)],[HDIlim(2,2)-mean(nontarget(:,2)) HDIlim(4,2)-mean(nontarget(:,4))])
else
errorbar([1 2 3 4],[mean(nontarget(:,1)) mean(nontarget(:,3)) mean(nontarget(:,5)) mean(nontarget(:,7))],[mean(nontarget(:,1))-HDIlim(1,1) mean(nontarget(:,3))-HDIlim(3,1) mean(nontarget(:,5))-HDIlim(5,1) mean(nontarget(:,7))-HDIlim(7,1)],[HDIlim(1,2)-mean(nontarget(:,1)) HDIlim(3,2)-mean(nontarget(:,3)) HDIlim(5,2)-mean(nontarget(:,5)) HDIlim(7,2)-mean(nontarget(:,7))])
errorbar([1 2 3 4],[mean(nontarget(:,2)) mean(nontarget(:,4)) mean(nontarget(:,6)) mean(nontarget(:,8))],[mean(nontarget(:,2))-HDIlim(2,1) mean(nontarget(:,4))-HDIlim(4,1) mean(nontarget(:,6))-HDIlim(6,1) mean(nontarget(:,8))-HDIlim(8,1)],[HDIlim(2,2)-mean(nontarget(:,2)) HDIlim(4,2)-mean(nontarget(:,4)) HDIlim(6,2)-mean(nontarget(:,6)) HDIlim(8,2)-mean(nontarget(:,8))])
end
hold off
ylim([0 1]);
xlim([0.8 (NC(EXPT)+0.2)]);
xticks([1:NC(EXPT)]);
legend(Legends);
if EXPT~=3
    xticklabels(Labels(EXPT,:));
end
title('Non-Target Response');

for i=1:2*NC(EXPT)
  HDIlim(i,:)=HDIofMCMC(uniform(:,i),0.95);
end

subplot(2,2,3);
hold on
if EXPT~=3
errorbar([1 2],[mean(uniform(:,1)) mean(uniform(:,3))],[mean(uniform(:,1))-HDIlim(1,1) mean(uniform(:,3))-HDIlim(3,1)],[HDIlim(1,2)-mean(uniform(:,1)) HDIlim(3,2)-mean(uniform(:,3))])
errorbar([1 2],[mean(uniform(:,2)) mean(uniform(:,4))],[mean(uniform(:,2))-HDIlim(2,1) mean(uniform(:,4))-HDIlim(4,1)],[HDIlim(2,2)-mean(uniform(:,2)) HDIlim(4,2)-mean(uniform(:,4))])
else
errorbar([1 2 3 4],[mean(uniform(:,1)) mean(uniform(:,3)) mean(uniform(:,5)) mean(uniform(:,7))],[mean(uniform(:,1))-HDIlim(1,1) mean(uniform(:,3))-HDIlim(3,1) mean(uniform(:,5))-HDIlim(5,1) mean(uniform(:,7))-HDIlim(7,1)],[HDIlim(1,2)-mean(uniform(:,1)) HDIlim(3,2)-mean(uniform(:,3)) HDIlim(5,2)-mean(uniform(:,5)) HDIlim(7,2)-mean(uniform(:,7))])
errorbar([1 2 3 4],[mean(uniform(:,2)) mean(uniform(:,4)) mean(uniform(:,6)) mean(uniform(:,8))],[mean(uniform(:,2))-HDIlim(2,1) mean(uniform(:,4))-HDIlim(4,1) mean(uniform(:,6))-HDIlim(6,1) mean(uniform(:,8))-HDIlim(8,1)],[HDIlim(2,2)-mean(uniform(:,2)) HDIlim(4,2)-mean(uniform(:,4)) HDIlim(6,2)-mean(uniform(:,6)) HDIlim(8,2)-mean(uniform(:,8))])
end
hold off
ylim([0 1]);
xlim([0.8 (NC(EXPT)+0.2)]);
xticks([1:NC(EXPT)]);
legend(Legends);
if EXPT~=3
    xticklabels(Labels(EXPT,:));
end
title('Uniform Response');

for i=1:2*NC(EXPT)
  HDIlim(i,:)=HDIofMCMC(precision(:,i),0.95);
end

subplot(2,2,4);
hold on
if EXPT~=3
errorbar([1 2],[mean(precision(:,1)) mean(precision(:,3))],[mean(precision(:,1))-HDIlim(1,1) mean(precision(:,3))-HDIlim(3,1)],[HDIlim(1,2)-mean(precision(:,1)) HDIlim(3,2)-mean(precision(:,3))])
errorbar([1 2],[mean(precision(:,2)) mean(precision(:,4))],[mean(precision(:,2))-HDIlim(2,1) mean(precision(:,4))-HDIlim(4,1)],[HDIlim(2,2)-mean(precision(:,2)) HDIlim(4,2)-mean(precision(:,4))])
else
errorbar([1 2 3 4],[mean(precision(:,1)) mean(precision(:,3)) mean(precision(:,5)) mean(precision(:,7))],[mean(precision(:,1))-HDIlim(1,1) mean(precision(:,3))-HDIlim(3,1) mean(precision(:,5))-HDIlim(5,1) mean(precision(:,7))-HDIlim(7,1)],[HDIlim(1,2)-mean(precision(:,1)) HDIlim(3,2)-mean(precision(:,3)) HDIlim(5,2)-mean(precision(:,5)) HDIlim(7,2)-mean(precision(:,7))])
errorbar([1 2 3 4],[mean(precision(:,2)) mean(precision(:,4)) mean(precision(:,6)) mean(precision(:,8))],[mean(precision(:,2))-HDIlim(2,1) mean(precision(:,4))-HDIlim(4,1) mean(precision(:,6))-HDIlim(6,1) mean(precision(:,8))-HDIlim(8,1)],[HDIlim(2,2)-mean(precision(:,2)) HDIlim(4,2)-mean(precision(:,4)) HDIlim(6,2)-mean(precision(:,6)) HDIlim(8,2)-mean(precision(:,8))])
end
hold off
ylim([0 maxprecision(EXPT)]);
xlim([0.8 (NC(EXPT)+0.2)]);
xticks([1:NC(EXPT)]);
legend(Legends);
if EXPT~=3
    xticklabels(Labels(EXPT,:));
end
title('Precision');


