function [] = line_plot(time,SBP,DBP,PulseR)
 MAP = (SBP + 2*DBP )/3;
 DP = PulseR .* SBP;
 
 figure('name','Analysis: Original data');    
    h = plot(time,SBP,time,DBP,time,PulseR);
    hold on;
        xlabel('x-axis');
        ylabel('y-axis');
        linspace(0,16);
        grid on
        grid minor;
         set(h(1),'linewidth',1.5);
         set(h(2),'linewidth',1.5);
         set(h(3),'linewidth',1.5);
        yticks(0:5:140);
        ylim([0,140]);
        xticks(0:1:16);
        xlim([min(time) max(time)]);
        x = time;
        ystart = DBP;
        ystop = SBP;
        tx = [x.';x.';nan(1,length(x))];
        ty = [ystart.';ystop.';nan(1,length(x))];
        plot(tx(:),ty(:),'b');
        legend(h(1:3),'Average SBP','Average DBP','Average PulseR');
    hold off;
    
    
 figure('name','Analysis: MAP'); 
    h2 = plot(time, MAP,'b');
    hold on;
        xlabel('time');
        ylabel('MAP');
        linspace(0,16);
        grid on
        grid minor;
            set(h2,'linewidth',1.5)
        yticks(0:5:140);
        ylim([0,140]);
        xticks(0:1:16);
        xlim([min(time) max(time)]);
        legend('MAP');
     hold off;
    
     
    
 figure('name','Analysis: Double Product');
    h2 = plot(time, DP,'b');
    hold on;
        xlabel('time');
        ylabel('Double Product');
        grid on
        grid minor;
            set(h2,'linewidth',1.5)
        yticks(0:500:max(DP));
        ylim([0 max(DP)]);
        xticks(0:1:16);
        xlim([min(time) max(time)]);
        legend('Double Product');
     hold off;
 
