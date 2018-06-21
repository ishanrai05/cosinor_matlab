function [CI_Amp_min1, CI_Amp_max1, CI_phi_min1, CI_phi_max1,CI_Amp_min2, CI_Amp_max2, CI_phi_min2, CI_phi_max2,CI_Amp_min3, CI_Amp_max3, CI_phi_min3, CI_phi_max3] = CIcalc(X,T,Z,beta1,gamma1,n,sigma1,Amp1,phi1,beta2,gamma2,sigma2,Amp2,phi2,beta3,gamma3, sigma3,Amp3,phi3,alpha)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Confidence Interval Calculations
%
% Description: 
%   SubFuction of 'cosinor.m'. Finds individual confidence intervals for 
%   Amplitude and Acrophase and plots a polar plot representation of the
%   cosinor fit.
%
%   Follows cosinor analysis of a time series as outlined by
%   Nelson et al. "Methods for Cosinor-Rhythmometry" Chronobiologica.
%   1979. Please consult reference.
%
% Parent Function:
%   'cosinor.m'
%
% Example: Run Parent Function
%   Define time series: 
%       y = [102,96.8,97,92.5,95,93,99.4,99.8,105.5];
%       t = [97,130,167.5,187.5,218,247.5,285,315,337.5]/360;
%   Define cycle length and alpha:
%       w = 2*pi;
%       alpha = .05;
%   Run Code:
%       cosinor(t,y,w,alpha)
%
% Record of revisions:
%     Date           Programmmer        Description of change
%     =====          ===========        ======================
%     6/10/08        Casey Cox          Original Code
%     6/24/08        Casey Cox          Revisions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Find beta and gamma confidence region
F_distr = finv(1-alpha,2,n-3);

    A = X;
    B = 2*T;
    C = Z;
    D1 = -2*X*beta1 - 2*T*gamma1;
    E1 = -2*T*beta1 - 2*Z*gamma1;
    F1 = X*beta1^2 + 2*T*beta1*gamma1 + Z*gamma1^2 - (2/n)*sigma1^2*F_distr;

    g_max1 = -(2*A*E1 - D1*B)/(4*A*C - B^2);
    
    D2 = -2*X*beta2 - 2*T*gamma2;
    E2 = -2*T*beta2 - 2*Z*gamma2;
    F2 = X*beta2^2 + 2*T*beta2*gamma2 + Z*gamma2^2 - (2/n)*sigma2^2*F_distr;

    g_max2 = -(2*A*E2 - D2*B)/(4*A*C - B^2);
    
    D3 = -2*X*beta3 - 2*T*gamma3;
    E3 = -2*T*beta3 - 2*Z*gamma3;
    F3 = X*beta3^2 + 2*T*beta3*gamma3 + Z*gamma3^2 - (2/n)*sigma3^2*F_distr;

    g_max3 = -(2*A*E3 - D3*B)/(4*A*C - B^2);

gamma_s1 = [g_max1-Amp1*2:Amp1/1000:g_max1+Amp1*2];
beta_s11 = (-(B.*gamma_s1 + D1) + sqrt((B.*gamma_s1 + D1).^2 - 4*A*(C.*gamma_s1.^2 + E1.*gamma_s1 + F1)))/(2*A);
beta_s21 = (-(B.*gamma_s1 + D1) - sqrt((B.*gamma_s1 + D1).^2 - 4*A*(C.*gamma_s1.^2 + E1.*gamma_s1 + F1)))/(2*A);

gamma_s2 = [g_max2-Amp2*2:Amp2/1000:g_max2+Amp2*2];
beta_s12 = (-(B.*gamma_s2 + D2) + sqrt((B.*gamma_s2 + D2).^2 - 4*A*(C.*gamma_s2.^2 + E2.*gamma_s2 + F2)))/(2*A);
beta_s22 = (-(B.*gamma_s2 + D2) - sqrt((B.*gamma_s2 + D2).^2 - 4*A*(C.*gamma_s2.^2 + E2.*gamma_s2 + F2)))/(2*A);


gamma_s3 = [g_max3-Amp3*2:Amp3/1000:g_max3+Amp3*2];
beta_s13 = (-(B.*gamma_s3 + D3) + sqrt((B.*gamma_s3 + D3).^2 - 4*A*(C.*gamma_s3.^2 + E3.*gamma_s3 + F3)))/(2*A);
beta_s23 = (-(B.*gamma_s3 + D3) - sqrt((B.*gamma_s3 + D3).^2 - 4*A*(C.*gamma_s3.^2 + E3.*gamma_s3 + F3)))/(2*A);


%Isolate ellipse region
IND1 = find(real(beta_s11) ~= real(beta_s21));
gamma_s1 = gamma_s1(IND1); beta_s11 = beta_s11(IND1); beta_s21 = beta_s21(IND1);

IND2 = find(real(beta_s12) ~= real(beta_s22));
gamma_s2 = gamma_s2(IND2); beta_s12 = beta_s12(IND2); beta_s22 = beta_s22(IND2);

IND3 = find(real(beta_s13) ~= real(beta_s23));
gamma_s3 = gamma_s3(IND3); beta_s13 = beta_s13(IND3); beta_s23 = beta_s23(IND3);

%Determine if confidence region overlaps the pole.
if (range(gamma_s1) >= max(gamma_s1)) && ((range(beta_s11) >= max(beta_s11)) || (range(beta_s21) >= max(beta_s21)))
    disp('!! Confidence region overlaps the pole. Confidence limits for Amplitude and Acrophase cannot be determined !!');disp(' ');
    
    CI_Amp_max1 = [0];
    CI_Amp_min1 = [0];
    CI_phi_max1 = [0];
    CI_phi_min1 = [0];
else
    %Confidence Intervals for Amplitude
    CI_Amp_max1 = max(max([sqrt(beta_s11.^2 + gamma_s1.^2); sqrt(beta_s21.^2 + gamma_s1.^2)],[],2));
    CI_Amp_min1 = min(min([sqrt(beta_s11.^2 + gamma_s1.^2); sqrt(beta_s21.^2 + gamma_s1.^2)],[],2));
    
    %Confidence Intervals for Acrophase
    theta1 = cat(2,[atan(abs(gamma_s1./beta_s11))], [atan(abs(gamma_s1./beta_s21))]);
        a1 = sign(cat(2,[beta_s11],[beta_s21]));
        b1 = sign(cat(2,[gamma_s1],[gamma_s1]))*3;
        c1 = a1 + b1;
        for ii = 1:length(c1);
            if (c1(ii) == 4 || c1(ii) == 3)
                CIphi1(ii) = -theta1(ii);
                c1(ii) = 1;
            elseif (c1(ii) == 2 || c1(ii) == -1) 
                CIphi1(ii) = -pi + theta1(ii);
                c1(ii) = 2;
            elseif (c1(ii) == -4 || c1(ii) == -3)
                CIphi1(ii) = -pi - theta1(ii);
                c1(ii) = 3;
            elseif (c1(ii) == -2 || c1(ii) == 1)
                CIphi1(ii) = -2*pi + theta1(ii);
                c1(ii) = 4;
            end
        end
    if max(c1) - min(c1) == 3   
        CI_phi_max1 = min(CIphi1(c == 1));
        CI_phi_min1 = max(CIphi1(c == 4));
    else
        CI_phi_max1 = max(CIphi1);
        CI_phi_min1 = min(CIphi1);
    end
end



if (range(gamma_s2) >= max(gamma_s2)) && ((range(beta_s12) >= max(beta_s12)) || (range(beta_s22) >= max(beta_s22)))
    disp('!! Confidence region overlaps the pole. Confidence limits for Amplitude and Acrophase cannot be determined !!');disp(' ');
    
    CI_Amp_max2 = [0];
    CI_Amp_min2 = [0];
    CI_phi_max2 = [0];
    CI_phi_min2 = [0];
else
    %Confidence Intervals for Amplitude
    CI_Amp_max2 = max(max([sqrt(beta_s12.^2 + gamma_s2.^2); sqrt(beta_s22.^2 + gamma_s2.^2)],[],2));
    CI_Amp_min2 = min(min([sqrt(beta_s12.^2 + gamma_s2.^2); sqrt(beta_s22.^2 + gamma_s2.^2)],[],2));
    
    %Confidence Intervals for Acrophase
    theta2 = cat(2,[atan(abs(gamma_s2./beta_s12))], [atan(abs(gamma_s2./beta_s22))]);
        a2 = sign(cat(2,[beta_s12],[beta_s22]));
        b2 = sign(cat(2,[gamma_s2],[gamma_s2]))*3;
        c2 = a2 + b2;
        for ii = 1:length(c2);
            if (c2(ii) == 4 || c2(ii) == 3)
                CIphi2(ii) = -theta2(ii);
                c2(ii) = 1;
            elseif (c2(ii) == 2 || c2(ii) == -1) 
                CIphi2(ii) = -pi + theta2(ii);
                c2(ii) = 2;
            elseif (c(ii) == -4 || c(ii) == -3)
                CIphi(ii) = -pi - theta(ii);
                c(ii) = 3;
            elseif (c2(ii) == -2 || c2(ii) == 1)
                CIphi2(ii) = -2*pi + theta2(ii);
                c2(ii) = 4;
            end
        end
    if max(c2) - min(c2) == 3   
        CI_phi_max2 = min(CIphi2(c2 == 1));
        CI_phi_min2 = max(CIphi2(c2 == 4));
    else
        CI_phi_max2 = max(CIphi2);
        CI_phi_min2 = min(CIphi2);
    end
end


if (range(gamma_s3) >= max(gamma_s3)) && ((range(beta_s13) >= max(beta_s13)) || (range(beta_s23) >= max(beta_s23)))
    disp('!! Confidence region overlaps the pole. Confidence limits for Amplitude and Acrophase cannot be determined !!');disp(' ');
    
    CI_Amp_max3 = [0];
    CI_Amp_min3 = [0];
    CI_phi_max3 = [0];
    CI_phi_min3 = [0];
else
    %Confidence Intervals for Amplitude
    CI_Amp_max3 = max(max([sqrt(beta_s13.^2 + gamma_s3.^2); sqrt(beta_s23.^2 + gamma_s3.^2)],[],2));
    CI_Amp_min3 = min(min([sqrt(beta_s13.^2 + gamma_s3.^2); sqrt(beta_s23.^2 + gamma_s3.^2)],[],2));
    
    %Confidence Intervals for Acrophase
    theta3 = cat(2,[atan(abs(gamma_s3./beta_s13))], [atan(abs(gamma_s3./beta_s23))]);
        a3 = sign(cat(2,[beta_s13],[beta_s23]));
        b3 = sign(cat(2,[gamma_s3],[gamma_s3]))*3;
        c3 = a3 + b3;
        for ii = 1:length(c3);
            if (c3(ii) == 4 || c3(ii) == 3)
                CIphi3(ii) = -theta3(ii);
                c3(ii) = 1;
            elseif (c3(ii) == 2 || c3(ii) == -1) 
                CIphi3(ii) = -pi + theta3(ii);
                c3(ii) = 2;
            elseif (c3(ii) == -4 || c3(ii) == -3)
                CIphi3(ii) = -pi - theta3(ii);
                c3(ii) = 3;
            elseif (c3(ii) == -2 || c3(ii) == 1)
                CIphi3(ii) = -2*pi + theta3(ii);
                c3(ii) = 4;
            end
        end
    if max(c3) - min(c3) == 3   
        CI_phi_max3 = min(CIphi3(c3 == 1));
        CI_phi_min3 = max(CIphi3(c3 == 4));
    else
        CI_phi_max3 = max(CIphi3);
        CI_phi_min3 = min(CIphi3);
    end
end

%Polar Representation of Cosinor analysis 
figure('name','Rhythm Parameter Estimates with Joint Confidence Region', 'position', [245 357 643 600]);
plot(gamma_s1,beta_s11,'linewidth', 2); hold on;
plot(gamma_s2,beta_s12,'linewidth', 2);
plot(gamma_s3,beta_s13,'linewidth', 2);
plot(gamma_s1,beta_s21,'linewidth', 2);
plot(gamma_s2,beta_s22,'linewidth', 2);
plot( gamma_s3,beta_s23,'linewidth', 2);
    line([0 gamma1], [0 beta1], 'color','k','linewidth', 2)
    line([0 gamma2], [0 beta2], 'color','k','linewidth', 2)
    line([0 gamma3], [0 beta3], 'color','k','linewidth', 2)
    line([gamma1 -2*Amp1*sin(phi1)], [beta1 2*Amp1*cos(phi1)], 'color','k','linewidth', 2, 'linestyle',':')
    line([gamma2 -2*Amp2*sin(phi2)], [beta2 2*Amp2*cos(phi2)], 'color','k','linewidth', 2, 'linestyle',':')
    line([gamma3 -2*Amp3*sin(phi3)], [beta3 2*Amp3*cos(phi3)], 'color','k','linewidth', 2, 'linestyle',':')
    xlabel('\gamma')
    ylabel('\beta')
    ylim([-Amp1*2.5 Amp1*2.5])
    xlim([-Amp1*2.5 Amp1*2.5])
    line([0 0], [-Amp1*2 Amp1*2], 'color','k','linestyle', '--')
    line([0 0], [-Amp2*2 Amp3*2], 'color','k','linestyle', '--')
    line([0 0], [-Amp3*2 Amp3*2], 'color','k','linestyle', '--')
    line([-Amp1*2 Amp1*2], [0 0], 'color','k','linestyle', '--')
    line([-Amp2*2 Amp2*2], [0 0], 'color','k','linestyle', '--')
    line([-Amp3*2 Amp3*2], [0 0], 'color','k','linestyle', '--')
    
    %Clock and Labels
    theta_clock = (0:pi/60:2*pi)';
    clock1 = ([Amp1*2*cos(theta_clock) Amp1*2*sin(theta_clock)]);
    clock2 = ([Amp2*2*cos(theta_clock) Amp2*2*sin(theta_clock)]);
    clock3 = ([Amp3*2*cos(theta_clock) Amp3*2*sin(theta_clock)]);
    plot(clock1(:,1),clock1(:,2), 'k', 'linewidth', 1.5)
    plot(clock2(:,1),clock2(:,2), 'k', 'linewidth', 1.5)
    plot(clock3(:,1),clock3(:,2), 'k', 'linewidth', 1.5)
    plot(clock1(:,1)*1.2,clock1(:,2)*1.2, 'k', 'linewidth', 1.5)
    plot(clock2(:,1)*1.2,clock2(:,2)*1.2, 'k', 'linewidth', 1.5)
    plot(clock3(:,1)*1.2,clock3(:,2)*1.2, 'k', 'linewidth', 1.5)
        theta_clock = (0:pi/4:2*pi-pi/4)';
        clock_labels_xy = ([Amp1*2.2*cos(theta_clock) Amp1*2.2*sin(theta_clock)]);
        clock_labels = {'6:00'; '3:00'; '0:00'; '21:00'; '18:00'; '15:00'; '12:00'; '9:00'};
        for jj=1:length(clock_labels)
            text(clock_labels_xy(jj,1), clock_labels_xy(jj,2), clock_labels(jj), 'horizontalalignment', 'center');
        end