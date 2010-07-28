function subjects = processCSVfiles(params)

% function to make plots from CSV files generated by EM Clustering
% Module. The input params passes the required parameters to the function:
%  params.fileDir             location where CSV files are stored
%  params.clusterIds          cluster Ids for which plots need to be generated

% The following optinal parameters can also be specified:
%  params.NAval               the integer value by which the missing data
%                             is represented (default -1)
%  params.sm                  smoothing factor which is applied to the
%                             plotted curves (default 5)

% Mahnaz Maddah, 2009.

params.NAval = -1;   % what represents missing data
params.sm = 1;       % smoothing factor
population = params.population;
%%%%%
fileDir = params.fileDir;
clusterIds = params.clusterIds;
NAval = params.NAval;
sm = params.sm;
subjects = params.subjects;

if ~strcmp('/',fileDir(end))
    fileDir = [fileDir,'/'];
end

for k= 1:length(clusterIds)

   % figure('Name',sprintf('Diffusion parameters versus normalized arc length for cluster %d',k),'NumberTitle','off')
    [FA] = load([fileDir, sprintf('cluster%d_FA.csv',clusterIds(k))]);           %NxS
    [MD] = load([fileDir, sprintf('cluster%d_MD.csv',clusterIds(k))]);           %NxS
    [parDiff] = load([fileDir, sprintf('cluster%d_ParDiff.csv',clusterIds(k))]); %NxS
    [preDiff] = load([fileDir, sprintf('cluster%d_PerDiff.csv',clusterIds(k))]); %NxS
    
    if population
        fid = fopen([fileDir, sprintf('cluster%d_subjectNames.txt', clusterIds(k))]);
        trajectories_subjectName = textscan(fid, '%s');
        trajectories_subjectName =trajectories_subjectName{1};
        
        for s=1:length(subjects)
            %return cell Ids of trajectories_subjectNames that are equal to subjects{s}.name
            cellIds = find( strcmp([trajectories_subjectName],[params.subjectsDirectory subjects(s).name]));
            
            if (~isempty(cellIds))

                subjects(s).data = cell(4,1);
            
                subjects(s).data{1} = FA(cellIds,:);
                subjects(s).data{2} = MD(cellIds,:);
                subjects(s).data{3} = parDiff(cellIds,:);
                subjects(s).data{4} = preDiff(cellIds,:);
                [subjects(s).data_means, subjects(s).data_stds, subjects(s).w] = getStatsScalarMeasures(subjects(s).data,NAval);
                plotScalarMeasures(subjects(s).data_means, subjects(s).data_stds, sm, subjects(s).group);
            else
                %disp('subject has no contributing trajectories');
                subjects(s).data =[];
            end
        end
    else    
    
    data =cell(4);
    data{1} = FA;
    data{2} = MD;
    data{3} = parDiff;
    data{4} = preDiff;
    
    [data_means, data_stds] = getStatsScalarMeasures(data,NAval);
    plotScalarMeasures(data_means, data_stds, sm);
   
    end
end
end

function c = smooth(varargin)
%Smooths data by moving average else
y = varargin{1};
y = y(:);   % make it a column

span = [];
if nargin > 1
    span = varargin{2};
end
if isempty(span)
    span = 5;
end

span = floor(span);
n = length(y);
span = min(span,n);
w = span-1+mod(span,2); % make it odd
if w==1,
    c = y;
    return;
end

c = filter(ones(w,1)/w,1,y);
cbegin = cumsum(y(1:w-2));
cbegin = cbegin(1:2:end)./(1:2:(w-2))';
cend = cumsum(y(n:-1:n-w+3));
cend = cend(end:-2:1)./(w-2:-2:1)';
c = [cbegin;c(w:end);cend];
end

function [data_means, data_stds, w] = getStatsScalarMeasures(alldata,NAval)

for d=1:4 %over features:
    data = alldata{d};
    w = zeros(1,size(data,2));
    for j=1:size(data,2)
        nonzero = find(data(:,j)~= NAval);
        w(j) = length(nonzero);
        if w(j)>0
            data_mean(j) = mean(data(nonzero,j));
            data_std(j) = std(data(nonzero,j));
        else
            data_mean(j)=NaN;
            data_std(j)=NaN;
        end
        
       
    end
    data_means(d,:) = data_mean;
    data_stds(d,:) = data_std;
end
end

function plotScalarMeasures(data_means, data_stds, sm, group)
% This function makes a single figure with four subplots of diffusion
% scalar measures, 'FA', 'MD','Parallel Diffusion' and 'Perpendicular
% Diffusion', stored in the fisrt input argument of this function.
% Features are smoothed to the extent specified by the second input.
features = [{'FA'},{'MD'},{'Parallel Diffusion'}, {'Perpendicular Diffusion'}];
hold on
col ='b';
if exist('group')
    if group ==1
        col = 'r';
    end
end
  
for d=1:4
    
    data_mean = data_means(d,:);
    data_std = data_stds(d,:);
    
    data_mean(find(isnan(data_mean)))=0;
    data_std(find(isnan(data_std)))=0;
    
    N = length(data_mean);
    subplot(2,2,d)
    plot(smooth(data_mean,sm),'.','MarkerFaceColor',col, 'MarkerEdgeColor',col);%'LineWidth',2);
    hold on
   % plot(smooth(data_mean+data_std,sm),'--b','LineWidth',2);
   % plot(smooth(data_mean-data_std,sm),'--b','LineWidth',2);
    set(gca,'XTick',0:N/2:N)
    set(gca,'XTickLabel',{'0', '0.5','1'});
    ylabel(features(d))
    xlim([0 N]);

end
end