%% 0:   Create vector with all the subject paths

subject = dir('Data/EEG_fMRI_2017/Raw_data');
% removes hidden files from list
subject(1:3) = [];

temp = extractfield(subject, 'folder');
stats.path = temp(:);

temp = extractfield(subject, 'name');
stats.raw_name = temp(:);

names = convertCharsToStrings(extractfield(subject, 'name'));
stats.ID = names(:);

stats.path = strcat(stats.path,'/' + stats.ID)



%%  1:  Load physiological variables (heart rate and respiration) and global signal (GS) from HCP data
HR_data = {};
HR_avgs = {};
GS_data = {}
resp_data = {}
paths = stats.path;
nScans = length(paths)

r_PRF_sc_data = []

paths = stats.path;
for i = 1:7
    a = strcat(stats.path,'/Phys_sum.mat');
    b = strcat(stats.path,'/TissueBasedRegressors.mat');
    load(a(i))
    load(b(i))



    sc = 140;     % choosing a scan (sc) from 1-164, 41 patients who have 4 scans each

    GS=zscore(WB.MA(:)); HR=HRV(:); resp=zscore(resp(:)); % rows is time, column is scan number
    Ts_10 = 0.1;                                                       % Sampling period in seconds
    time_10 = 0:Ts_10:(length(HR)-1)*Ts_10;    % getting the time series, start from 0 and go up by ts until the scaled version of the last value
    % timeMR = time_10(ind_BOLD_10);            % indices on where you got an image from the MRI; time_10 is sampling for 10Hz



    HR_data = [HR_data; HR];
    HR_avgs = [HR_avgs; mean(HR)];
    GS_data = [GS_data; GS];
    resp_data = [resp_data; resp];
    close all




%% 2: Estimate PRF_sc parameters ***







    
%derive TIMEindices for mri; Time_10, time_MR
ind_BOLD_10 = zeros(length(timeMR),1);
NV = length(timeMR);
for j = 1:NV
    [M,I] = min(abs(time_10 - timeMR(j)));
    
    ind_BOLD_10(j) = I;
end


% ***response resp data; why 10*1.5
resp_s = smooth(resp_10,10*1.5);
% RF is derivative, add a zero and the rest and square everything
RF=diff(resp_s); RF=[0;RF(:)]; RF = RF.^2;


% genetic algorithm  generate with options
ga_opts = gaoptimset('TolFun',1e-10,'StallGenLimit',20,'Generations',30,'Display','iter','UseParallel',1);   % Display: iter
options = optimoptions('fmincon','Display','off','Algorithm','interior-point',...
    'UseParallel',true,'MaxIterations',60,'MaxFunctionEvaluations',3000,'OptimalityTolerance',1e-8,'PlotFcn','optimplotfval');    % 'PlotFcn','optimplotfval'

x0 = [  1    2.5   5.6    0.9    1.9   2.9   12.5    0.5 ];  
ub = x0+3;
lb = x0-3; lb(find(lb<0))=0;

% HR:heart rate; RF:respiratory flow, ind_BOLD_10:indices for
% downsampling,GS: global signal
% *why put @(p)
h = @(P) PRF_sc_optimize_parameters(P,Ts_10,HR,RF,ind_BOLD_10,GS);
%     x0 = ga(h,length(ub),[],[],[],[],lb,ub,[],[],ga_opts);
x_opt = fmincon(h,x0,[],[],[],[],lb,ub,[],options);

[obj_function,CRF_sc,RRF_sc,HR_conv,RF_conv,r_PRF_sc,yPred] = h(x_opt);

fprintf(' ----------------------------------------------- \n')
fprintf('Correlation b/w GS and PRF output \n')
fprintf('CRF (HR): %3.1f%%  \n',r_PRF_sc(2)*100)
fprintf('RRF (RF): %3.1f%%  \n',r_PRF_sc(3)*100)
fprintf('CRF & RRF (HR & RF): %3.1f%%  \n',r_PRF_sc(1)*100)



r_PRF_sc_data =  [r_PRF_sc_data; r_PRF_sc(:)'];
end

stats.r_PRF_sc_data = r_PRF_sc_data

disp(mean(stats.r_PRF_sc_data(1)))


stats.HR_data = HR_data
stats.HR_avg = HR_avgs
stats.GS = GS_data
stats.resp = resp_data
disp(mean(stats.r_PRF_sc_data(:,1)))