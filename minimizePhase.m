function [phase,dphase,allStats] = minimizePhase(x,levs1,sem1,levs2,sem2,nMonte,doplot,figNum)

% x1 = x coordinate of levs1( in nt) 

delta = 0.01;
phaseTest= -1:delta:1;

allStats.nMonte = nMonte;
allStats.delta = delta;
allStats.phaseScore = cell(1,nMonte);
allStats.spline1 = cell(1,nMonte);
allStats.spline2 = cell(1,nMonte);
allStats.phase = zeros(1,nMonte);

for ii = 1:nMonte

phaseScore = 0*phaseTest;
xSpline = min(x):delta:max(x); 

levs1M = levs1 + randn(1,length(levs1)).*sem1;
levs2M = levs2 + randn(1,length(levs2)).*sem2;

allStats.spline1{ii} = spline(x,levs1M,xSpline);
allStats.spline2{ii} = spline(x,levs2M,xSpline);

for jj = 1:length(phaseScore)
    s2temp = spline(x,levs2M,xSpline-phaseTest(jj));
    score = -(s2temp - allStats.spline1{ii}).^2;
    phaseScore(jj) = sum(score);
end
[~,maxInd] = max(phaseScore);
allStats.phaseScore{ii} = phaseScore;
allStats.phase(ii) = phaseTest(maxInd);
    
end

if doplot
    figure(figNum); clf ; 
        subplot(3,1,1);hold on
            plot(xSpline,allStats.spline1{end},'color',jcolor(12),'linewidth',2)
            plot(xSpline,allStats.spline2{end},'color',jcolor(14),'linewidth',2)
            plot(xSpline+mean(allStats.phase),allStats.spline2{end},...
                'color','k','linewidth',2,'linestyle','--')
        subplot(3,1,2);hold on
            h = histogram(allStats.phase,'normalization','pdf','binlimits',[min(phaseTest) max(phaseTest)],'binwidth',0.05);
            plot(phaseTest,normpdf(phaseTest,mean(allStats.phase),std(allStats.phase)),'r','linewidth',2)
        subplot(3,1,3);hold on
            plot(phaseTest,phaseScore,'linewidth',2)
end

phase = mean(allStats.phase);
dphase = std(allStats.phase);

end

