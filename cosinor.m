function [] = cosinor(t,SBP,DBP,PulseR,w,alpha)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COSINOR 	[]=cosinor(t,y,w,alpha)
%
% Description: 
%   Cosinor analysis uses the least squares method to fit a sine wave to a
%   time series. Cosinor analysis is often used in the analysis
%   of biologic time series that demonstrate predictible rhythms. This
%   method can be used with an unequally spaced time series.
%
%   Follows cosinor analysis of a time series as outlined by
%   Nelson et al. "Methods for Cosinor-Rhythmometry" Chronobiologica.
%   1979. Please consult reference.
%
% Input:
%   t - time series
%   y - value of series at time t
%   w - cycle length, defined by user based on prior knowledge of time
%       series
%   alpha - type I error used for cofidence interval calculations. Usually 
%       set to be 0.05 which corresponds with 95% cofidence intervals
%
% Define Variables:
%   M - Mesor, the average cylce value
%   Amp - Amplitude, half the distance between peaks of the fitted
%       waveform
%   phi - Acrophase, time point in the cycle of highest amplitude (in
%       radians)
%   RSS - Residual Sum of Squares, a measure of the deviation of the
%       cosinor fit from the original waveform
%
% Subfunctions:
%   'CIcalc.m'
%
% Example:
%   Define time series: 
%       y = [102,96.8,97,92.5,95,93,99.4,99.8,105.5];
%       t = [97,130,167.5,187.5,218,247.5,285,315,337.5]/360;
%   Define cycle length and alpha:
%       w = 2*pi;
%       alpha = .05;
%   Run Code:
%       cosinor(t,y,w,alpha)

% Record of revisions:
%     Date           Programmmer        Description of change
%     =====          ===========        ======================
%     5/16/08        Casey Cox          Original Code
%     6/24/08        Casey Cox          Revisions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin ~= 6
    error('Incorrect number of inputs');
end
if length(t) < 4
    error('There must be atleast four time measurements')
end

%% I. Parameter Estimation

n = length(t);

% Substituition
x = cos(w.*t);
z = sin(w.*t);

% Set up and solve the normal equations simultaneously
NE1 = [  n        sum(x)       sum(z)     sum(SBP);
      sum(x)   sum(x.^2)    sum(x.*z)    sum(x.*SBP);
      sum(z)   sum(x.*z)    sum(z.^2)    sum(z.*SBP);];

RNE1 = rref(NE1);
NE2 = [  n        sum(x)       sum(z)     sum(DBP);
      sum(x)   sum(x.^2)    sum(x.*z)    sum(x.*DBP);
      sum(z)   sum(x.*z)    sum(z.^2)    sum(z.*DBP);];

RNE2 = rref(NE2);
NE3 = [  n        sum(x)       sum(z)     sum(PulseR);
      sum(x)   sum(x.^2)    sum(x.*z)    sum(x.*PulseR);
      sum(z)   sum(x.*z)    sum(z.^2)    sum(z.*PulseR);];

RNE3 = rref(NE3);
M1 = RNE1(1,4); beta1 = RNE1(2,4); gamma1 = RNE1(3,4);
M2 = RNE2(1,4); beta2 = RNE2(2,4); gamma2 = RNE2(3,4);
M3 = RNE3(1,4); beta3 = RNE3(2,4); gamma3 = RNE3(3,4);

%Calculate amplitude and acrophase from beta and gamma
Amp1 = sqrt(beta1^2 + gamma1^2);
theta1 = atan(abs(gamma1/beta1));

Amp2 = sqrt(beta2^2 + gamma2^2);
theta2 = atan(abs(gamma2/beta2));

Amp3 = sqrt(beta3^2 + gamma3^2);
theta3 = atan(abs(gamma3/beta3));

    % Calculate acrophase (phi) and convert from radians to degrees
    a1 = sign(beta1);
    b1 = sign(gamma1);
    if (a1 == 1 || a1 == 0) && b1 == 1
        phi1 = -theta1;
    elseif a1 == -1 && (b1 == 1 || b1 == 0) 
        phi1 = -pi + theta1;
    elseif (a1 == -1 || a1 == 0) && b1 == -1
        phi1 = -pi - theta1;
    elseif a1 == 1 && (b1 == -1 || b1 == 0)
        phi1 = -2*pi + theta1;
    end
    
    a2 = sign(beta2);
    b2 = sign(gamma2);
    if (a2 == 1 || a2 == 0) && b2 == 1
        phi2 = -theta2;
    elseif a2 == -1 && (b2 == 1 || b2 == 0) 
        phi2 = -pi + theta2;
    elseif (a2 == -1 || a2 == 0) && b2 == -1
        phi2 = -pi - theta2;
    elseif a2 == 1 && (b2 == -1 || b2 == 0)
        phi2 = -2*pi + theta2;
    end
    
    a3 = sign(beta3);
    b3 = sign(gamma3);
    if (a3 == 1 || a3 == 0) && b3 == 1
        phi3 = -theta3;
    elseif a3 == -1 && (b3 == 1 || b3 == 0) 
        phi3 = -pi + theta3;
    elseif (a3 == -1 || a3 == 0) && b3 == -1
        phi3 = -pi - theta3;
    elseif a3 == 1 && (b3 == -1 || b3 == 0)
        phi3 = -2*pi + theta3;
    end

% Display results
disp('Parameters:'); disp('---------------');
fprintf(1,'Mesor = %g \nAmplitude = %g \nAcrophase = %g \n\n',M1,Amp1,phi1);
disp('Parameters:'); disp('---------------');
fprintf(1,'Mesor = %g \nAmplitude = %g \nAcrophase = %g \n\n',M2,Amp2,phi2);

disp('Parameters:'); disp('---------------');
fprintf(1,'Mesor = %g \nAmplitude = %g \nAcrophase = %g \n\n',M3,Amp3,phi3);


%Plot orginal data and cosine fit
f1 = M1 + Amp1*cos(w.*t+phi1);
f2 = M2 + Amp2*cos(w.*t+phi2);
f3 = M3 + Amp3*cos(w.*t+phi3);

    figure('name','Cosinor Analysis: Original data and fitted function');
    plot(t,SBP,t,DBP,t,PulseR); hold on;
        xlabel('x-axis');
        ylabel('y-axis');
        grid on;
        grid minor;
    h = plot(t,f1,'r',t,f2,'b',t,f3,'g');
         set(h(1),'linewidth',2);
         set(h(2),'linewidth',2);
         set(h(3),'linewidth',2);
         grid on;
         grid minor;
        legend('Original', 'Cosinor');
        xlim([min(t) max(t)]);
 
        

%% II. Confidence Limtes for Single Cosinor

%Residual sum of errors
RSS1 = sum((SBP - (M1 + beta1.*x + gamma1.*z)).^2);
RSS2 = sum((DBP - (M2 + beta2.*x + gamma2.*z)).^2);
RSS3 = sum((PulseR - (M3 + beta3.*x + gamma3.*z)).^2);

%Residual varience estimation
sigma1 = sqrt(RSS1/(n-3));
sigma2 = sqrt(RSS2/(n-3));
sigma3 = sqrt(RSS3/(n-3));

%Find confidence interval for mesor
    X = 1/n * sum((x - mean(x)).^2);
    Z = 1/n * sum((z - mean(z)).^2);
    T = 1/n * sum((x - mean(x)).*(z - mean(z)));

%Confidence interval for the mesor
CI_M1 = tinv(1-alpha/2,n-3)*sigma1^2*sqrt(((sum(x.^2))*(sum(z.^2)) - (sum(x.*z))^2)/(n^3*(X*Z - T^2))); %#ok<NASGU>
CI_M2 = tinv(1-alpha/2,n-3)*sigma2^2*sqrt(((sum(x.^2))*(sum(z.^2)) - (sum(x.*z))^2)/(n^3*(X*Z - T^2))); %#ok<NASGU>
CI_M3 = tinv(1-alpha/2,n-3)*sigma3^2*sqrt(((sum(x.^2))*(sum(z.^2)) - (sum(x.*z))^2)/(n^3*(X*Z - T^2))); %#ok<NASGU>

%Find confidence intervals for the amplitude and acrophase
[CI_Amp_min1, CI_Amp_max1, CI_phi_min1, CI_phi_max1,CI_Amp_min2, CI_Amp_max2, CI_phi_min2, CI_phi_max2,CI_Amp_min3, CI_Amp_max3, CI_phi_min3, CI_phi_max3] = CIcalc(X,T,Z,beta1,gamma1,n,sigma1,Amp1,phi1,beta2,gamma2,sigma2,Amp2,phi2,beta3,gamma3,sigma3,Amp3,phi3,alpha); %#ok<NASGU,NASGU>



%% III. Zero-amplitude test
p_3a1 = fpdf((n*(X*beta1^2 + 2*T*beta1*gamma1 + Z*gamma1^2)/(2*sigma1^2)),2,n-3);
p_3a2 = fpdf((n*(X*beta2^2 + 2*T*beta2*gamma2 + Z*gamma2^2)/(2*sigma2^2)),2,n-3);
p_3a3 = fpdf((n*(X*beta3^2 + 2*T*beta3*gamma3 + Z*gamma3^2)/(2*sigma3^2)),2,n-3);
fprintf(1,'Zero Amplitude Test \n')
fprintf(1,'------------------------------------------------------\n')
fprintf(1,'Amplitude        0.95 Confidence Limits        P Value\n')
fprintf(1,'---------        ----------------------        -------\n')
fprintf(1,' %.2f               (%.2f to %.2f)             %g\n\n',Amp1,CI_Amp_min1,CI_Amp_max1,p_3a1)
fprintf(1,' %.2f               (%.2f to %.2f)             %g\n\n',Amp2,CI_Amp_min2,CI_Amp_max2,p_3a2)
fprintf(1,' %.2f               (%.2f to %.2f)             %g\n\n',Amp3,CI_Amp_min3,CI_Amp_max3,p_3a3)
