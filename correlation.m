function [] = correlation(SBP,DBP,PulseR)

 figure('name','SBP?DBP');    
 
 hold on;
        plotregression(DBP,SBP,'Regression');
        yticks(0:10:140);
        ylim([0,140]);
        xticks(0:10:140);
        xlim([0,140]);
    hold off;

 figure('name','Analysis: Original data');    
    scatter (PulseR,SBP,0.1);
    hold on;
        plotregression(PulseR,SBP,'Regression');
        yticks(0:10:140);
        ylim([0,140]);
        xticks(0:10:140);
        xlim([0,140]);
    hold off;
   