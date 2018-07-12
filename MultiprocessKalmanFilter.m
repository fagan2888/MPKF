function [y_smooth, p_j, p_pj, p_ppj] = MultiprocessKalmanFilter(y_in, t_in, p0, Ki, m0, C0)

% This function is a first attempt to construct the Extended Kalman Filter
% used by the Lattec/Delaval Herd Navigator to smooth data and assign
% trends to the data
% The input for this filter is a time series of values y_in e.g. milk P4
% taken at timepoints t_in
%       y_in = raw time series
%       t_in = time points at which raw time series is taken
%       p0   = 4 parameter vector of initial probabilities to be in a
%              certain state
%       K    = 4x3 matrix containing the variance-correcting parameters
%              K(i,:) = [Kv Kg Kd] for each state i
%              Kv = measurement error variance multiplier if state i 
%                   1*c if j = 1, 2, 3 and large if i = 4 (outlier)
%              Kg = if the level changes, Kg is large  and µt is adapted
%              Kd = if the slope changes, Kd is large and betat is adapted
%              always proportional to c = measurement error constant 
% The output consists of a time of corrected milk P4 values y_out and the
% state characteristics, containing posterior probabilities, one step back
% probabilities and two steps back probabilities. These outputs might be 
% supplemented with additional outputs later
%       y_out = smoothed time series
%       M1 = series of posterior probabilities for state 1 + one and two
%            step back probabilities
%       M2 = series of posterior probabilities for state 2 + one and two
%            step back probabilities
%       M3 = series of posterior probabilities for state 3 + one and two
%            step back probabilities
%       M4 = series of posterior probabilities for state 4 + one and two
%            step back probabilities

% Assumption 1: the probability of the previous observation to be in state
%               i does not affect the probability of the current 
%               observation to be in state j. <-> Patent: "Model j is
%               assumed to be selected with probability P(Mt(j)|Dt-1) =
%               p0(j), independently of the past Dt-1, j=1,...,4 (fixed
%               model selection probabilities)
% Assumption 2: A priori it is assumed that teta0 ~ N(m0,C0) and that all
%               of the parameters are known (i.e. p0, K, m0, C0)
% Assumption 3: The actual sampling moment is not taken into account, only
%               the rank/order. So, we use the n-th measurement in stead of
%               the actual time point. This might be untrue, and later
%               evaluations should point out whether the distance between
%               measurements can be taken into account to estimate beta_t
%                   Possibly, time interval t can be taken into account by
%                   using r(t)/ r(t-1) * beta(t-1) with r(t) length of time
%                   interval between t and t-1 
%                   According to patent, not implemented, because Gt is
%                   taken as G = [1 1; 0 1]; otherwise, Gt would be equal
%                   to Gt = [1 r(t)/r(t-1); 0 r(t)/r(t-1)];
% Assumption 4: c is the measurement error, proportional to the absolute
%               value of the time series inputed. We assume this c² to be
%               known and equal to 1, having no effect. In a later stage,
%               we might need c if it turns out that the variance depends
%               on the level. In principle, one would assume that this is
%               the case for P4, but having the convergence of values as 28
%               ng/mL, we don't know this is implemented?
% Assumption 5: The guess of p0 will be strongly dependent on the sampling
%               scheme; if the samples are taken 'smart', another sampling
%               p0 will apply then when samples are taken over the whole
%               lactation. In papers: p0 is determined based on an
%               historical dataset of 100 cows taken over 6 months. No idea
%               whether this is adjusted per farm or over time. We will
%               start assuming a very broad range of values for p0, and
%               test this on the 'normal' bio-model determined sampling.
%               if later it turns out that the parameters are not valid
%               anymore for the increased sampling, we will estimate new
%               parameters here. This complicates the comparisons 


%% Time-series characterization & clean-up
% Simulate data for script testing
    % clear all
    % close all
    % load('Sim_data.mat')

    % y_in = [4 + 3.*randn(9,1)' , -20+3*randn(15,1)' , -560+21.6.*(25:50)+3*randn(26,1)'];
    % y_in(40) = 330;
    % save('Sim_data.mat','y_in')

	% t_in = 1:50;
    % figure; 
    % plot(t_in,y_in,'LineWidth',2); hold on; plot(t_in,y_in,'.','MarkerSize',12)

    % m0 = [4 0]';                % initial concentration (P4) = 1 ng/mL
    % C0 = [20 0;0 10];       % best guess for relatively great uncertainty about initial mean/intercept and trend

    % assumptions on initial parameters p0 for each model state j
    % taken from Korsgaard and Lovendahl (2002) might be totally wrong
    % p0 = [0.94 0.02 0.02 0.02];         % prior chance y is in state 1, 2, 3 of 4 
                                    % [SS; change level; change slope; outlier];
                            
    % assumptions on variance multipliers K = [Kv Kg Kd] for each state j
    % K = [1 0 0; 1 20 0; 1 0 10; 50 0 0];
    % K = [9 0 0; 9 20 0; 9 0 10; 50 0 0];


% check whether length of input data are equal
L = length(y_in);
    if length(t_in) ~= L
            error('Unequal length of y_in and t_in')
    end

% check for empty values
ind = find(isnan(y_in) == 0);

Y_in = y_in(ind);
T_in = t_in(ind);                   % reset to start from zero

L = length(Y_in);                   % length of input array for recursive updating

% R = diff(T_in);                   % unequally spaced observations taken into account





%% Linear Growth model
% we assume that the time-series can be described with a linear growth
% model at each time t

F = [1 0];              % F is the observation matrix of the obs equation
G = [1 1; 0 1];         % G is the evolution matrix of the system equation


% construct variance multiplier from the input matrix K corresponding to each state j
% Ki = [Kv Kd Kg Ko]
K = [Ki(1); Ki(1); Ki(1); Ki(1)+Ki(4)];

Kw.J1 = [0 0 ;0 0];         % state 1 = steady state
Kw.J2 = [Ki(2) 0 ; 0 0];    % state 2 = change of level
Kw.J3 = [Ki(2)+Ki(3) Ki(3) ; Ki(3) Ki(3 )];    % state 3 = change of slope
Kw.J4 = [0 0 ;0 0];    % state 4 = an outlier


% assume c² = known, equal to 1;

c = 1;                  % c = constant representing normal measurement error
                        % typically, this is the CV of the analyser used,
                        % which is a multiplicative measurement error.
                        % Might be a problem when using HN data


%% start modelling at t = 1, 
% using the assumption that system at time 0 is normally distributed with
% mean 0 and variance c²*C0
% goal: smooth time series Y_in, given the probabilities that previous
% measurments were in previous states and the prior information available
% to be calculated per step:
%   p_j, p_pj, p_ppj vectors of length L, L-1 and L-2 containing the probabilities of the
%                  system to be in state 1-4 at time t, time t-1 and time
%                  t-2
% the smoothed time-series is only calculated at time t, but not updated
% later when the probabilities of time t-1 are updated. However, these count
% when calculating variance-covariance matrix to estimate system equation

% time = 1
t = 1;

% Step 1: construct the backbones of the needed variables - forecast of y
y_out(1,:) = zeros(1,L);        % output vector contains same number of values as the input vector
y_out(2,:) = T_in;              % time points at which y_out is calculated

% forecast y1 - this is forecast but not smoothed value.
y_out(1,t) = F*G*m0;

% for t = 1; we assume that at time 0, that the initial assumption of the distribution of
% teta0 ~ N(m0,c²C0), then the likelihood that y1 at time 1 belongs to model j is
% given by p(y1 | M1(j)) is obtained from 
% (y1|M(j),M(i),D0) ~ N(F*G*m0(i)  ,  c^2*(Kv(j)*mu^2 + F*(G*C0(i)G'+Kw(j))*F'))
% by filling in m0 and C0   % eq. 3.12
% % % p(1:5,t) =0;
py(1:4,t) = 0;
for j = 1:4
    mod = sprintf('J%d',j);
        % the prob that y_in is a measurement of state j = likelihood
    py(j,t) = normpdf(Y_in(t), F*G*m0 , sqrt(c^2*(K(j,1)+F*(G*C0*G'+Kw.(mod))*F'))); % likelyhood
        % the prob that y_in is a measurement of state j * prior
        % probability that state j is valid
% % %     p(j,t) = py(j,t)*p0(j); % starting value of p(j=1)
end


% calculate p_ij ---> if we assume steady state, this should be p(1,t)?? 
% Calculate likelihood that y is a measurement of model i when at t-1
% y(t-1) was a measurement of model j
% rows : j ; cols : i
p_ij = zeros(4,4);      % 3.10
for j = 1:4                     % model j at time 1 = ROWS
    for i = 1:4                 % model i at time 0 = COLS
        p_ij(j,i) = py(j,t)*p0(i)*p0(j);
    end
end

% normalization
norm = 1/sum(sum(p_ij));

p_ij = p_ij*norm;

p_j(:,t) = sum(p_ij,2);         % probability to be in state j at time t
    
% p_pj(:,t) = sum(p_ij)';       % probability to be in state j at time t-1


% calculate posterior distribution of teta1 (teta1|M1(j),D1) ~N(m1(j),C1(j)
% to this end, mixture probabilities are taken into consideration
% step 1: calculate eq. 3.8 C(i,j)
    % in the final step: use m(j-1) for µ!!!! 
C_ij = zeros(8,8);
for j = 1:4             % ROW = current model
    for i = 1:4         % COL = former model
        mod = sprintf('J%d',j);
%         C_ij(2*j-1:2*j,2*i-1:2*i) = (F'*(K(j,1)*m0(1)^2)^(-1)*F + (G*C0*G'+Kw.(mod))^(-1))^(-1);     % j = 1:4
        C_ij(2*j-1:2*j,2*i-1:2*i) = (F'*(K(j,1))^(-1)*F + (G*C0*G'+Kw.(mod))^(-1))^(-1);     % j = 1:4

    end
end


% step 2: calculate m_ij at time t=1 (eq.3.9) -- final: use m(j-1) --- 3.16!!
m_ij = zeros(8,4);
for j = 1:4             % ROW = current model
    for i = 1:4         % COL = former model
%         m_ij(2*j-1:2*j,i) = C_ij(2*j-1:2*j,2*i-1:2*i)*F'*(K(j,1)*m0(1))^(-1)*(Y_in(t)-F*G*m0)+G*m0;     % j = 1:4
            m_ij(2*j-1:2*j,i) = C_ij(2*j-1:2*j,2*i-1:2*i)*F'*(K(j,1))^(-1)*(Y_in(t)-F*G*m0)+G*m0;     % j = 1:4
    end
end

% calculate m_j at time t=1 (eq.3.14)
m_j = zeros(8,L);
for j = 1:4             % ROW = current model
    for i = 1:4         % COL = former model
        m_j(2*j-1:2*j,t) = m_j(2*j-1:2*j,t) + p_ij(j,i)/p_j(j,t)*m_ij(2*j-1:2*j,i);     % m_j is [level(t); slope(t)] for model j at time t; eq.3.14
%         m_j(2*j-1:2*j,t) = m_j(2*j-1:2*j,t) + (p_ij2(i,j)*p(5,t))*m_ij(2*j-1:2*j,i);     % m_j is [level(t); slope(t)] for model j at time t; eq.3.14
    end
end

% calculate C_j at time t=1 (eq.3.15)
C_j = zeros(8,2*L);
for j = 1:4             % ROW = current model
    for i = 1:4         % COL = former model
        C_j(2*j-1:2*j,2*t-1:2*t) = C_j(2*j-1:2*j,2*t-1:2*t) + p_ij(j,i)/p_j(j,t)*(C_ij(2*j-1:2*j,2*i-1:2*i) + (m_ij(2*j-1:2*j,i)-m_j(2*j-1:2*j,t))*(m_ij(2*j-1:2*j,i)-m_j(2*j-1:2*j,t))');    % m_j is [level(t); slope(t)] for model j at time t; eq.3.14
%         C_j(2*j-1:2*j,2*t-1:2*t) = C_j(2*j-1:2*j,2*t-1:2*t) + p_ij2(i,j)*p(5,t)*(C_ij(2*j-1:2*j,2*i-1:2*i) + (m_ij(2*j-1:2*j,i)-m_j(2*j-1:2*j,t))*(m_ij(2*j-1:2*j,i)-m_j(2*j-1:2*j,t))');    % m_j is [level(t); slope(t)] for model j at time t; eq.3.14
    end
end


% calculate m at time t to fill in, in step 2, in the variance of syst eq acc to eq.
% 3.16, first component -- not sure of this step!! % reflects not the right
% number, why would this be so small at this time??
m = zeros(2,L);

y_smooth(1,t) = 0;
for j = 1:4
    y_smooth(1,t) = y_smooth(1,t) + p_j(j,t)*F*m_j(2*j-1:2*j,t);     % smoothed ?
end
y_smooth(2,1:L) = Y_in;


%% recursive updating procedure, giving the value y_out each time t and the posterior and one and two step back probabilities that model j is valid at time t
% figure; 
% subplot(1,2,1); hold on
% plot(y_out(2,:),Y_in,'.-')
% % plot(y_out(2,t),y_out(1,t),'ro')
% plot(y_out(2,t),y_smooth(1,t),'ro')
% legend({'Raw','Corrected'},'AutoUpdate','off')
% 
% subplot(3,2,2); hold on;
% plot(t,p_j(1,t) ,'.','Color',[0 0.5 0],'Linewidth',1,'MarkerSize',10)
% plot(t,p_j(2,t) ,'s','Color',[1 0   1],'Linewidth',2,'MarkerSize',5)
% plot(t,p_j(3,t) ,'^','Color',[0 0   1],'Linewidth',2,'MarkerSize',5)
% plot(t,p_j(4,t) ,'o','Color',[1 0.5 0],'Linewidth',2,'MarkerSize',5)
% legend({'Steady state','Level change','Slope change','Outlier'},'AutoUpdate','off')


for t = 2:L
    
    % calculate m(t-1) eq. 3.16 and forecast of y_out 
    for j = 1:4
        m(1:2,t-1) = m(1:2,t-1) + p_j(j,t-1)*m_j(2*j-1:2*j,t-1);  %eq. 3.16:
        y_out(1,t) = F*G*m(1:2,t-1);     %eq. 3.16: forecast ?
    end
    
%     subplot(1,2,1); hold on
%     plot(y_out(2,t),y_out(1,t),'ro')
%     plot(t,m(1,t-1),'ro')    

    % calculate forecast distribution of yt conditional on M(j) gives the
    % probability of p(yt|Mt(j),Mt-1(i),Dt-1), needed to calculate p_ij
    % likelihood
    py_ij = zeros(4,4);
    for j = 1:4             % ROW = current model
        mod = sprintf('J%d',j);
        for i = 1:4         % COL = former model
%             py_ij(j,i) = normpdf(Y_in(t), F*G*m_j(2*i-1:2*i,t-1),sqrt(c^2*(K(j,1)*m(1,t-1)^2+F*(G*C_j(2*i-1:2*i,2*(t-1)-1:2*(t-1))*G'+Kw.(mod))*F')));
            py_ij(j,i) = normpdf(Y_in(t), F*G*m_j(2*i-1:2*i,t-1),sqrt(c^2*(K(j,1)+F*(G*C_j(2*i-1:2*i,2*(t-1)-1:2*(t-1))*G'+Kw.(mod))*F')));
        end
    end
    
    % store for two step back calculations
    p_pij = p_ij;
    
    % calculate p_ij at time t (eq.3.10)
    % here, we take pt-1,(i) = p0(i)
    for j = 1:4             % ROW = current model
        for i = 1:4         % COL = former model
            p_ij(j,i) = py_ij(j,i)*p_j(i,t-1)*p0(j);
        end
    end
    
    norm = 1/sum(sum(p_ij));
    p_ij = p_ij*norm;
    p_j(:,t) = sum(p_ij,2);
    p_pj(:,t-1) = sum(p_ij)';


%     subplot(3,2,2); hold on;
%         plot(t,p_j(1,t) ,'.','Color',[0 0.5 0],'Linewidth',1,'MarkerSize',10)
%         plot(t,p_j(2,t) ,'s','Color',[1 0   1],'Linewidth',2,'MarkerSize',5)
%         plot(t,p_j(3,t) ,'^','Color',[0 0   1],'Linewidth',2,'MarkerSize',5)
%         plot(t,p_j(4,t) ,'o','Color',[1 0.5 0],'Linewidth',2,'MarkerSize',5)


%     subplot(3,2,4); hold on;
%         plot(t-1,p_pj(1,t-1) ,'.','Color',[0 0.5 0],'Linewidth',1,'MarkerSize',10)
%         plot(t-1,p_pj(2,t-1) ,'s','Color',[1 0   1],'Linewidth',2,'MarkerSize',5)
%         plot(t-1,p_pj(3,t-1) ,'^','Color',[0 0   1],'Linewidth',2,'MarkerSize',5)
%         plot(t-1,p_pj(4,t-1) ,'o','Color',[1 0.5 0],'Linewidth',2,'MarkerSize',5)
%         if t == 2; legend({'Steady state','Level change','Slope change','Outlier'},'AutoUpdate','off'); end
    
    % calculate posterior distribution of teta at time t via 3.8, 3.9, 3.14, 3.15
        % calculate C_ij (eq. 3.8)
    for j = 1:4             % ROW = current model
        for i = 1:4         % COL = former model
            mod = sprintf('J%d',j);
%             C_ij(2*j-1:2*j,2*i-1:2*i) = (F'*(K(j,1)*m(1,t-1)^2)^(-1)*F + (G*C_j(2*i-1:2*i,2*(t-1)-1:2*(t-1))*G'+Kw.(mod))^(-1))^(-1);     % j = 1:4
            C_ij(2*j-1:2*j,2*i-1:2*i) = (F'*(K(j,1))^(-1)*F + (G*C_j(2*i-1:2*i,2*(t-1)-1:2*(t-1))*G'+Kw.(mod))^(-1))^(-1);     % j = 1:4
            
        end
    end
    
    % calculate m_ij (eq. 3.9)
    for j = 1:4             % ROW = current model
        for i = 1:4         % COL = former model
%             m_ij(2*j-1:2*j,i) = C_ij(2*j-1:2*j,2*i-1:2*i)*F'*(K(j,1)*m(1,t-1)^2)^(-1)*(Y_in(t)-F*G*m_j(2*i-1:2*i,t-1))+G*m_j(2*i-1:2*i,t-1);     % j = 1:4
             m_ij(2*j-1:2*j,i) = C_ij(2*j-1:2*j,2*i-1:2*i)*F'*(K(j,1))^(-1)*(Y_in(t)-F*G*m_j(2*i-1:2*i,t-1))+G*m_j(2*i-1:2*i,t-1);     % j = 1:4
        end
    end

    % calculate m_j at time t (eq. 3.14)
    m_j(:,t) = 0;
    for j = 1:4             % ROW = current model
        for i = 1:4         % COL = former model
            m_j(2*j-1:2*j,t) = m_j(2*j-1:2*j,t) + p_ij(j,i)/p_j(j,t)*m_ij(2*j-1:2*j,i);     % m_j is [level(t); slope(t)] for model j at time t; eq.3.14
%             m_j(2*j-1:2*j,t) = m_j(2*j-1:2*j,t) + p_ij2(j,i)*p(j,t)*m_ij(2*j-1:2*j,i);     % m_j is [level(t); slope(t)] for model j at time t; eq.3.14
        end
    end
   
    % calculate C_j at time t (eq. 3.15)
    C_j(:,2*t-1:2*t) = 0;
    for j = 1:4             % ROW = current model
        for i = 1:4         % COL = former model
            C_j(2*j-1:2*j,2*t-1:2*t) = C_j(2*j-1:2*j,2*t-1:2*t) + p_ij(j,i)/p_j(j,t)*(C_ij(2*j-1:2*j,2*i-1:2*i) + (m_ij(2*j-1:2*j,i)-m_j(2*j-1:2*j,t))*(m_ij(2*j-1:2*j,i)-m_j(2*j-1:2*j,t))');    % m_j is [level(t); slope(t)] for model j at time t; eq.3.14
        end
    end
    
    
    % calculate 1 step back probability of state j at time t, given info on t+1: p(M(j)|Dt,yt+1)
    y_smooth(1,t) = 0;
    for j = 1:4
        y_smooth(1,t) = y_smooth(1,t) + p_j(j,t)*F*m_j(2*j-1:2*j,t);     %eq. 3.16: forecast ?
    end
%     subplot(1,2,1)
%     plot(y_out(2,t),y_smooth(1,t),'ro')
    
    % calculate 2 step back probabilities of state j at time t given info on time t+2 and t+1 (p(M(j)|Dt,yt+2)
    if t > 2
        py_ki(4,4) = 0;
        for i = 1:4         % ROW = current model 
            mod = sprintf('J%d',i);
            for k = 1:4     % COL = one step back previous model
                py_ki(i,k) = normpdf(Y_in(t), F*G*m_j(2*k-1:2*k,t-1),sqrt(c^2*(K(i,1)+F*(G*C_j(2*k-1:2*k,2*(t-1)-1:2*(t-1))*G'+Kw.(mod))*F')));
            end
        end
        
        for i = 1:4         % two steps back previous model
            a(1:4,1:4,i)= py_ki*p0(i);
        end
        
        for i = 1:4
            a(1:4,1:4,i) = a(1:4,1:4,i)*p_pij;
        end
        norm = 1/sum(sum(sum(a)));
        a = a*norm;
        
        % sum over k
        for j = 1:4
            b(1:4,j) = sum(a(1:4,1:4,j))';
        end
        
        % sum over i
        p_ppj(:,t-2) = sum(b,2)';

%         subplot(3,2,6); hold on;
%             plot(t-2,p_ppj(1,t-2) ,'.','Color',[0 0.5 0],'Linewidth',1,'MarkerSize',10)
%             plot(t-2,p_ppj(2,t-2) ,'s','Color',[1 0   1],'Linewidth',2,'MarkerSize',5)
%             plot(t-2,p_ppj(3,t-2) ,'^','Color',[0 0   1],'Linewidth',2,'MarkerSize',5)
%             plot(t-2,p_ppj(4,t-2) ,'o','Color',[1 0.5 0],'Linewidth',2,'MarkerSize',5)
%             if t == 3; legend({'Steady state','Level change','Slope change','Outlier'},'AutoUpdate','off'); end
    end
end
y_smooth = y_smooth(1,:)';

% subplot(1,2,1)
%     plot(y_out(2,:),y_smooth(1,:),'r.-')
% title('Raw and smoothed values')
% xlabel('Time')
% ylabel('Measurements')
% 
% subplot(3,2,2); ylim([0 1])
%     plot(t_in,p_j(1,:) ,'.-','Color',[0 0.5 0],'Linewidth',1,'MarkerSize',10)
%     plot(t_in,p_j(2,:) ,'s-','Color',[1 0   1],'Linewidth',1,'MarkerSize',5)
%     plot(t_in,p_j(3,:) ,'^-','Color',[0 0   1],'Linewidth',1,'MarkerSize',5)
%     plot(t_in,p_j(4,:) ,'o-.','Color',[1 0.5 0],'Linewidth',1,'MarkerSize',5)
% title('Posterior probabilities')
% xlabel('Time')
% ylabel('Probability')
% 
% subplot(3,2,4); ylim([0 1])
%     plot(t_in(1:end-1),p_pj(1,:) ,'.-','Color',[0 0.5 0],'Linewidth',1,'MarkerSize',10)
%     plot(t_in(1:end-1),p_pj(2,:) ,'s-','Color',[1 0   1],'Linewidth',1,'MarkerSize',5)
%     plot(t_in(1:end-1),p_pj(3,:) ,'^-','Color',[0 0   1],'Linewidth',1,'MarkerSize',5)
%     plot(t_in(1:end-1),p_pj(4,:) ,'o-.','Color',[1 0.5 0],'Linewidth',1,'MarkerSize',5)
% title('One step back smoothed probabilities')
% xlabel('Time')
% ylabel('Probability')
% 
% subplot(3,2,6); ylim([0 1])
%     plot(t_in(1:end-2),p_ppj(1,:) ,'.-','Color',[0 0.5 0],'Linewidth',1,'MarkerSize',10)
%     plot(t_in(1:end-2),p_ppj(2,:) ,'s-','Color',[1 0   1],'Linewidth',1,'MarkerSize',5)
%     plot(t_in(1:end-2),p_ppj(3,:) ,'^-','Color',[0 0   1],'Linewidth',1,'MarkerSize',5)
%     plot(t_in(1:end-2),p_ppj(4,:) ,'o-.','Color',[1 0.5 0],'Linewidth',1,'MarkerSize',5)
% title('Two steps back smoothed probabilities')
% xlabel('Time')
% ylabel('Probability')

