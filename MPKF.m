function [y_smooth, p_j, p_pj, p_ppj] = MPKF(y_in, t_in, p0, Ki, m0, C0)
 
% This function contains the recursive updating procedure of a multiprocess, 
% ‘Extended’ Kalman Filter. To smooth progesterone time series, we suggest to use 
% p0 = [0.6966 0.20220.000022  0.1011] and Ki = [12 29 0.3 50] and have m0 dependent on
% the first measurement, C0 = [0.1 0; 0 0.1]
% The input for this filter is a time series of values y_in e.g. milk P4
% taken at time points t_in
%       y_in = raw time series
%       t_in = time points at which raw time series is taken
%       p0   = 4 parameter vector of initial probabilities to be in a
%              certain state
%       K    = 4x4 matrix containing the variance-correcting parameters
%              K(i,:) = [Kv Kg Kd Ko] for each state i
%              Kv = measurement error variance multiplier if state i 
%                   1*c if j = 1, 2, 3 and large if i = 4 (outlier)
%              Kg = if the level changes, Kg is large  and µt is adapted
%              Kd = if the slope changes, Kd is large and betat is adapted
%              always proportional to c = measurement error constant 
%		Ko = outlier variance multiplier
% The output consists of a time of corrected milk P4 values y_smooth and the
% state characteristics, containing posterior probabilities p_j, one step back
% probabilities p_pj and two steps back probabilities p_ppj.
%       y_smooth	= smoothed time series
%       p_j 		= series of posterior probabilities of each model j
%       p_pj 		= one step back probabilities of each model j
%       p_ppj 		= two steps back probabilities of each model j
%            
% Assumption 1: the probability of the previous observation to be in state
%               i does not affect the probability of the current 
%               observation to be in state j. <-> Patent: "Model j is
%               assumed to be selected with probability P(Mt(j)|Dt-1) =
%               p0(j), independently of the past Dt-1, j=1,...,4 (fixed
%               model selection probabilities)
% Assumption 2: A priori it is assumed that teta0 ~ N(m0,C0) and that all
%               of the parameters are known (i.e. p0, K, m0, C0)
% Assumption 3: The actual sampling moment is not taken into account, only
%               the rank/order. So, we use the n-th measurement instead of
%               the actual time point. 
%                   	Possibly, time interval t can be taken into account by
%                   	using r(t)/ r(t-1) * beta(t-1) with r(t) length of time
%                   	interval between t and t-1 
%                   	According to the patent (Friggens et al.2007), not implemented,
%	            	because Gt is taken as G = [1 1; 0 1]; otherwise, 
%			Gt would be equal to Gt = [1 r(t)/r(t-1); 0 r(t)/r(t-1)];
% Assumption 4: c, often introduced for proportial measurement errors, is not 
% 		 included absolute because the variance V = Kv is assumed constant
% Assumption 5: The guess of p0 will be strongly dependent on the sampling
%               scheme; if the samples are taken 'smart', another sampling
%               p0 will apply then when samples are taken over the whole
%               lactation. In papers: p0 is determined based on an
%               historical dataset of 100 cows taken over 6 months. No idea
%               whether this is adjusted per farm or over time. 
 
%% Time-series characterization & clean-up
% check length of input data are equal
L = length(y_in);
    if length(t_in) ~= L
        error('Unequal length of y_in and t_in')
    end
 
% check for empty values
Y_in = y_in(isnan(y_in) == 0);
T_in = t_in(isnan(y_in) == 0);      % reset to start from zero
 
L = length(Y_in);                   % length of input array for recursive updating
 

%% Linear Growth model
% we assume that the time-series can be described with a linear growth
% model at each time t
 
F = [1 0];              % F is the observation matrix of the observation equation
G = [1 1; 0 1];         % G is the evolution matrix of the system equation
 
% construct variance multiplier from the input matrix K corresponding to each state j
% Ki = [Kv Kd Kg Ko]
K = [Ki(1); Ki(1); Ki(1); Ki(1)+Ki(4)];
 	
Kw.J1 = [0 0 ;0 0];         		% state 1 = steady state
Kw.J2 = [Ki(2) 0 ; 0 0];    		% state 2 = change of level
Kw.J3 = [Ki(2)+Ki(3) Ki(3) ; Ki(3) Ki(3 )];    % state 3 = change of slope
Kw.J4 = [0 0 ;0 0];    		% state 4 = an outlier
 

%% Start modelling at t = 1, 
% goal: smooth time series Y_in, given the probabilities that previous
% measurements were in previous states and the prior information available
% to be calculated per step:
%   p_j, p_pj, p_ppj vectors of length L, L-1 and L-2 containing the probabilities of the
%                  system to be in state 1-4 at time t, time t-1 and time t-2    
% the smoothed time-series is only calculated at time t, but not updated
% later when the probabilities of time t-1 are updated. 
 
t = 1; 					% time = 1

% Step 1: construct the backbones of the needed variables - forecast of y
y_out(1,:) = zeros(1,L);        	% same number of values as the input vector
y_out(2,:) = T_in;              	% time points at which y_out is calculated
m = zeros(2,L); 			% m

% forecast y1 - this is forecast but not smoothed value.
y_out(1,t) = F*G*m0;
 
% for t = 1; we assume that at time 0, that the initial assumption of the distribution of
% teta0 ~ N(m0,C0), then the likelihood that y1 at time 1 belongs to model j is
% given by p(y1 | M1(j)) is obtained from 
% (y1|M(j),M(i),D0) ~ N(F*G*m0(i) , (Kv(j)+F*(G*C0(i)G'+Kw(j))*F'))

% by filling in m0 and C0   		% eq. S.25
py(1:4,t) = 0;				% py(j) contains ‘likelihood’ that y is observation of 						  state j
for j = 1:4
    mod = sprintf('J%d',j);		% model j
    py(j,t) = normpdf(Y_in(t), F*G*m0 , sqrt((K(j,1)+F*(G*C0*G'+Kw.(mod))*F'))); % likelihood
end
 
% calculate p_ij based on likelihood that y was a measurement of state i when at t-1
p_ij = zeros(4,4);      		% Eq. S.18
for j = 1:4                     	% model j at time 1 = ROWS
    for i = 1:4                 	% model i at time 0 = COLS
        p_ij(j,i) = py(j,t)*p0(i)*p0(j);
    end
end
 
% normalization step (Eq. S.19)
norm = 1/sum(sum(p_ij));
 
p_ij = p_ij*norm;			% Eq. S.18
p_j(:,t) = sum(p_ij,2);         	% Eq. S.20 - probability to be in state j at time t
p_pj(:,t) = sum(p_ij)';       	% Eq. S.29 - probability to be in state i at time t-1
 
 
% calculate posterior distribution of teta1 (teta1|M1(j),D1) ~N(m1(j),C1(j)
% to this end, mixture probabilities are taken into consideration
% step 1: calculate Eq. S.16 = C(i,j)
C_ij = zeros(8,8);
for j = 1:4             		% ROW = current model at time 1
    for i = 1:4         		% COL = former model at time 0 (prior)
        mod = sprintf('J%d',j);
        C_ij(2*j-1:2*j,2*i-1:2*i) = (F'*(K(j,1))^(-1)*F + (G*C0*G'+Kw.(mod))^(-1))^(-1);
    end
end
 
% step 2: calculate m_ij at time t=1 (eq.S.17) 
m_ij = zeros(8,4);
for j = 1:4             		% ROW = current model at time t
    for i = 1:4         		% COL = former model at time 0 (prior)
        m_ij(2*j-1:2*j,i) = C_ij(2*j-1:2*j,2*i-1:2*i)*F'*(K(j,1))^(-1)*(Y_in(t)-F*G*m0)+G*m0;     
    end
end
 
% calculate m_j at time t=1 (Eq. S.22)
m_j = zeros(8,L);
for j = 1:4             		% ROW = current model at time t
    for i = 1:4         		% COL = former model at time 0 (prior)
        m_j(2*j-1:2*j,t) = m_j(2*j-1:2*j,t) + p_ij(j,i)/p_j(j,t)*m_ij(2*j-1:2*j,i);     
        % m_j is [level(t); slope(t)] for model j at time t
    end
end
 
% calculate C_j at time t=1 (Eq. S.23)
C_j = zeros(8,2*L);
for j = 1:4             		% ROW = current model at time t
    for i = 1:4         		% COL = former model at time 0 (prior)        
        C_j(2*j-1:2*j,2*t-1:2*t) = C_j(2*j-1:2*j,2*t-1:2*t) + p_ij(j,i)/p_j(j,t)*(C_ij(2*j-1:2*j,2*i-1:2*i) + ...
                                   (m_ij(2*j-1:2*j,i)-m_j(2*j-1:2*j,t))*(m_ij(2*j-1:2*j,i)-m_j(2*j-1:2*j,t))');    		
    end
end
 
% calculate y_smooth (Eq. S.28)
y_smooth(1,t) = 0;
for j = 1:4
    y_smooth(1,t) = y_smooth(1,t) + p_j(j,t)*F*m_j(2*j-1:2*j,t);     
end
y_smooth(2,1:L) = Y_in;
 
 
%% Recursive updating procedure 
% giving the value y_smooth each time t and the posterior and one and two step back probabilities that model j is valid at time t

for t = 2:L
    
    % calculate m(t-1) and forecast of y_out 
    for j = 1:4
        m(1:2,t-1) = m(1:2,t-1) + p_j(j,t-1)*m_j(2*j-1:2*j,t-1);  
        y_out(1,t) = F*G*m(1:2,t-1);         
    end
    
    % calculate forecast distribution of yt conditional on M(j) given the
    % probability of p(yt|Mt(j),Mt-1(i),Dt-1), needed to calculate p_ij likelihood
    py_ij = zeros(4,4);
    for j = 1:4             		% ROW = current model at time t
        mod = sprintf('J%d',j);
        for i = 1:4         		% COL = former model at time t-1
            py_ij(j,i) = normpdf(Y_in(t), F*G*m_j(2*i-1:2*i,t-1),...
                         sqrt(c^2*(K(j,1)+F*(G*C_j(2*i-1:2*i,2*(t-1)-1:2*(t-1))*G'+Kw.(mod))*F')));
        end
    end
    
    % store for two step back calculations
    p_pij = p_ij;
    
    % calculate p_ij at time t (eq. S.18) with pt-1(i) = p0(i)
    for j = 1:4             		% ROW = current model at time t
        for i = 1:4         		% COL = former model at time t-1
            p_ij(j,i) = py_ij(j,i)*p_j(i,t-1)*p0(j);
        end
    end
    
    norm = 1/sum(sum(p_ij));		% Eq. S.19
    p_ij = p_ij*norm;			% Eq. S.18
    p_j(:,t) = sum(p_ij,2);		% Eq. S.20
    p_pj(:,t-1) = sum(p_ij)';		% Eq. S.29
 
     
    % calculate C_ij (Eq. S.16)
    for j = 1:4             		% ROW = current model at time t
        for i = 1:4         		% COL = former model at time t-1
            mod = sprintf('J%d',j);
            C_ij(2*j-1:2*j,2*i-1:2*i) = (F'*(K(j,1))^(-1)*F + ...
                                        (G*C_j(2*i-1:2*i,2*(t-1)-1:2*(t-1))*G'+Kw.(mod))^(-1))^(-1);    
        end
    end
    
    % calculate m_ij (Eq. S.17)
    for j = 1:4             		% ROW = current model at time t
        for i = 1:4         		% COL = former model at time t-1
             m_ij(2*j-1:2*j,i) = C_ij(2*j-1:2*j,2*i-1:2*i)*F'* ...
                                 (K(j,1))^(-1)*(Y_in(t)-F*G*m_j(2*i-1:2*i,t-1))+G*m_j(2*i-1:2*i,t-1);
        end
    end
 
    % calculate m_j at time t (Eq. S.22)
    m_j(:,t) = 0;
    for j = 1:4             		% ROW = current model at time t
        for i = 1:4         		% COL = former model at time t-1
            m_j(2*j-1:2*j,t) = m_j(2*j-1:2*j,t) + p_ij(j,i)/p_j(j,t)*m_ij(2*j-1:2*j,i);    
            % m_j is [level(t); slope(t)] for model j at time t
        end
    end
   
    % calculate C_j at time t (Eq. S.23)
    C_j(:,2*t-1:2*t) = 0;
    for j = 1:4             		% ROW = current model at time t
        for i = 1:4         		% COL = former model at time t-1
            C_j(2*j-1:2*j,2*t-1:2*t) = C_j(2*j-1:2*j,2*t-1:2*t) + ...
                                       p_ij(j,i)/p_j(j,t)*(C_ij(2*j-1:2*j,2*i-1:2*i) + ...
                                       (m_ij(2*j-1:2*j,i)-m_j(2*j-1:2*j,t))*(m_ij(2*j-1:2*j,i)-...
                                       m_j(2*j-1:2*j,t))');
        end
    end
    
    
    % calculate smoothed value (Eq. S.28)
    y_smooth(1,t) = 0;
    for j = 1:4
        y_smooth(1,t) = y_smooth(1,t) + p_j(j,t)*F*m_j(2*j-1:2*j,t);     
    end
    
    % calculate two step back probabilities of state j at time t given info on time t+2 and t+1
    if t > 2
        py_ki(4,4) = 0;		% likelihood
        for i = 1:4         		% ROW = current model at time t
            mod = sprintf('J%d',i);
            for k = 1:4     		% COL = previous model at time t-1
                py_ki(i,k) = normpdf(Y_in(t), F*G*m_j(2*k-1:2*k,t-1), sqrt(c^2*(K(i,1)+F*(G*C_j(2*k-1:2*k,2*(t-1)-1:2*(t-1))*G'+Kw.(mod))*F')));
            end
        end
        
        for i = 1:4         		% previous model at time t-2
            a(1:4,1:4,i)= py_ki*p0(i);
        end
        
        for i = 1:4			% normalization
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
     end
end

y_smooth = y_smooth(1,:)';
