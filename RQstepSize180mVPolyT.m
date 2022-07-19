

load('con3p')
load('con5p')

% visualize the consensus
close all
figure(1);subplot(1,2,1); ylabel('Current'); title('3p feeding'); xlabel('enzyme step number')
stairs(con3p.levels,'linewidth',1.25,'color','k')
subplot(1,2,2); ylabel('Current'); title('5p feeding'); xlabel('enzyme step number')
stairs(con5p.levels,'linewidth',1.25,'color','k')


% do the monte carlo & phase minimization
nMonte = 50;
warning('off','all')
% matlab doesn't like some of the NaN's here, but it doesn't really matter
levs.levs3pI = con3p.levels(1:2:35);
levs.stds3pI = con3p.sem(1:2:35);
levs.levs3pD = con3p.levels(2:2:36);
levs.stds3pD = con3p.sem(2:2:36);

levs.levs5pI = con5p.levels(2:2:20);
levs.stds5pI = con5p.sem(2:2:20);
levs.levs5pD = con5p.levels(3:2:21);
levs.stds5pD = con5p.sem(3:2:21);

levs.x3 = 1:length(levs.levs3pI);
[levs.phase3p,levs.dphase3p] = minimizePhase(levs.x3,levs.levs3pI,levs.stds3pI,...
    levs.levs3pD,levs.stds3pD,nMonte,true,1);

levs.x5 = 1:length(levs.levs5pI);
[levs.phase5p,levs.dphase5p] = minimizePhase(levs.x5,levs.levs5pI,levs.stds5pI,...
    levs.levs5pD,levs.stds5pD,nMonte,true,2);
warning('on','all')
%
figure(1); 
subplot(3,1,1)
title('Force assisting')
ylabel('Current (pA)')
ylim([10 60])
xlim([3 15])
set(gca,'xtick',3:15,'xticklabel',num2cell(0:12))
set(gca,'ytick',0:10:60)
xlabel('Position (nt)')
subplot(3,1,2);
ylabel('Probability (arbs)')
set(gca,'yticklabel','')
xlabel('Phase shift (nt)')
subplot(3,1,3);
ylabel('Score')
xlabel('Phase shift (nt)')
set(gca,'yscale','log')

figure(2); 
subplot(3,1,1)
title('Force opposing')
ylabel('Current (pA)')
xlabel('Position (nt)')
xlim([0 12])
ylim([10 60])
set(gca,'ytick',0:10:60,'xtick',0:12,'xticklabel',num2cell(0:12))
subplot(3,1,2)
ylabel('Probability  (arbs)')
set(gca,'yticklabel','')
xlabel('Phase shift (nt)')

subplot(3,1,3)
ylabel('Score')
xlabel('Phase shift (nt)')
set(gca,'yscale','log')

