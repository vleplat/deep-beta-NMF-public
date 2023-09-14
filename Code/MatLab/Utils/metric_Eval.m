function [SDR,SIR,SAR,perm] = metric_Eval(myFolder,extension,options)

% Add functions from Utils directory
addpath('./Utils');

% % Metrics computation 
x_se=[]; %source estimates (temporal)
filePattern = fullfile(myFolder, extension); % Change to whatever pattern you need.
se_theFiles = dir(filePattern);
for k = 1 : length(se_theFiles)
    baseFileName = se_theFiles(k).name;
    fullFileName = fullfile(myFolder, baseFileName);
    [x_se(k,:),options.fs]=audioread(fullFileName);
end

x_ts=[]; %true source (temporal)
myFolder = './Metric_eval/true_sources';
filePattern = fullfile(myFolder, extension); % Change to whatever pattern you need.
ts_theFiles = dir(filePattern);
if(length(se_theFiles)==length(ts_theFiles))
    for k = 1 : length(ts_theFiles)
        baseFileName = ts_theFiles(k).name;
        fullFileName = fullfile(myFolder, baseFileName);
        [x_ts(k,:),options.fs]=audioread(fullFileName);
    end
    [~,ln_se]=size(x_se);
    [~,ln_ts]=size(x_ts);
    if(ln_se<ln_ts)
        x_ts=x_ts(:,1:ln_se);
    elseif(ln_se>ln_ts)
        x_se=x_se(:,1:ln_ts);
    end
    [SDR,SIR,SAR,perm]=bss_eval_sources(x_se,x_ts);
else
    warning('number of files for source estimates and true source different');
end