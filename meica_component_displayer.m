function meica_component_displayer()

%select the text file that has the timecourses of interest
ctab = spm_select(Inf, '.*ctab.txt', 'Select ctab file');
ctab_loc = strfind(ctab,'ctab.txt');

ctab_loc =ctab_loc -1;
timecourses = strcat(ctab(1:ctab_loc),'mmix.1D');


fid = fopen(ctab);
tline = fgetl(fid);

while ischar(tline)
    %disp(tline)
    if strfind(tline, '#ACC')
        accp_list = tline;
    elseif strfind(tline, '#REJ')
        rej_list = tline;
    elseif strfind(tline, '#MID')
        mid_list = tline;
    elseif strfind(tline, '#IGN')
        ign_list = tline;
    elseif strfind(tline, '(VEx)');
        total_var = tline;
    end
    tline = fgetl(fid);
end
fclose(fid);

accp_list = accp_list(6:end);
rej_list = rej_list(6:end);
mid_list = mid_list(6:end);
ign_list = ign_list(6:end);
total_var = total_var(43:end);

accp_list = accp_list(1:(find(accp_list=='#')-1));
rej_list = rej_list(1:(find(rej_list=='#')-1));
mid_list = mid_list(1:(find(mid_list=='#')-1));
ign_list = ign_list(1:(find(ign_list=='#')-1));

accps = str2num(accp_list);
rejs = str2num(rej_list);
mids = str2num(mid_list);
igns = str2num(ign_list);

%Here we will turn these into useful matrices for comparison later.

%% Initialize variables.
filename = ctab;
delimiter = '\t';
startRow = 17;

%% Read columns of data as strings:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric strings to numbers.
% Replace non-numeric strings with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2,3,4,5,6]
    % Converts strings in the input cell array to numbers. Replaced non-numeric
    % strings with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(thousandsRegExp, ',', 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric strings to numbers.
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end


%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
imported_ctab = cell2mat(raw);
%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me R;

%timecourses = spm_select(Inf, '.*.1D', 'Select ctab file');

%For the imported_ctab, the columns are
%Comp#  Kappa   Rho   Variance   Normed_variance

savedir = spm_select(1,'dir','Select the MEICA output folder...');
cd(savedir);

ted_dir = strcat(savedir, '\TED');
cd(ted_dir); %Go into the output directory

all_betas = load_nii('betas_OC.nii');

all_betas = all_betas.img;

cd(savedir); %Go back out into the main directory. 

%%Making the motion plots, including framewise displacement
cfg.motionparam = 'motion.1D'; %output from MEICA, organized as: roll pitch yaw dS  dL  dP
cfg.prepro_suite = 'meica';
cfg.radius = 50;

figure;
grid on; grid minor;
[fwd,~]=bramila_framewiseDisplacement(cfg); %calculate FD using script

x_axis = size(fwd,1); %Get the number of timepoints

raw_motion = load(cfg.motionparam); %Load up motion for plotting.

%Subplots are used here to keep everything on the same screen.
% The y axes for the 6 motion estimates are set from the min and max of
% those values.
% The FD y axes is set to top out at 3 - as that is our voxel size.
% The purpose is two fold, anything larger than that is worrisome and
% This makes it easy to jump through all the subjects and compare
% quickly.
subplot(3,1,1); plot(raw_motion(:,4:6)); axis([0 x_axis min(min(raw_motion(:,4:6))) max(max(raw_motion(:,4:6)))]);
title('translation'); grid on;

rots = (raw_motion(:,1:3)); %These are already in degrees
subplot(3,1,2); plot(rots); axis([0 x_axis min(min(rots)) max(max(rots))]);
title('rotation'); grid on;

subplot(3,1,3); plot(fwd); title('Framewise Displacement'); axis([0 x_axis 0 3]);
grid on;

mkdir('component_plots');
savedir = [savedir, '/component_plots/'];

%plots are created here.
[~, titl, ~] = fileparts(cfg.motionparam);
print([savedir, titl], '-dpng');

timecourses_data = load(timecourses);


x_axis = size(timecourses_data,1);
figure;

%This is preparing for the brain slice display steps
num_betas = size(all_betas,4);

sag_slices = (size(all_betas,1));
cor_slices = (size(all_betas,2));
hor_slices = size(all_betas,3);

sag_cuts = floor(sag_slices/10);
cor_cuts = floor(cor_slices/10);
hor_cuts = floor(hor_slices/10);

color_table = [];

for i = 1:size(timecourses_data, 2);
    
    
    upper = max(max(max(all_betas(:,:,:,i))));
    lower = min(min(min(all_betas(:,:,:,i))));
    bounds = (0.3 * max([abs(upper), abs(lower)]));
    
    kappa = imported_ctab(i,2);
    rho = imported_ctab(i,3);
    variance_explained = imported_ctab(i,4);
    
    %plot(timecourses_data(:,i));
    %% Prints the graphs so that the outcome can be seen
    subplot(9,5,1:15)
    if any(accps == (i-1))
        plot(timecourses_data(:,i), 'Color',[0 .5 0]); %Green for BOLD like
        color_table = vertcat(color_table, [0 1 0]);
    elseif any(rejs == (i-1))
        plot(timecourses_data(:,i), 'r'); %Red for Rejected non BOLD
         color_table = vertcat(color_table, [1 0 0]);
    elseif any(mids == (i-1))
        plot(timecourses_data(:,i), 'm'); %Magenta for R2* weighted artifacts
         color_table = vertcat(color_table, [1 0 1]);
    elseif any(igns == (i-1))
        plot(timecourses_data(:,i), 'k'); %Black for Ignored components
         color_table = vertcat(color_table, [0 0 0]);
    end
    %%
    
    axis([0 x_axis min(timecourses_data(:,i)) max(timecourses_data(:,i))]);
      title(strcat('Component:', num2str(i), ', on ctab: ', num2str(i-1), ', kappa: ', num2str(kappa,3), ', rho: ', num2str(rho,3), ', variance: ', num2str(variance_explained,4))); grid on;
    label = strcat('Component_', num2str(i), '_on_ctab_', num2str(i-1));
    
    sag_img = [];
    for j = 1:sag_cuts:sag_slices
    sag_img = horzcat(sag_img, rot90(squeeze(all_betas(j,:,:, i))));
    end
    
    cor_img = [];
    for j = 1:cor_cuts:cor_slices
    cor_img = horzcat(cor_img, rot90(squeeze(all_betas(:,j,:, i))));
    end
    
    hor_img = [];
    for j = 1:hor_cuts:hor_slices
    hor_img = horzcat(hor_img, rot90(squeeze(all_betas(:,:,j, i))));
    end
    
    subplot(9,5,16:25)
    imshow(sag_img,[-bounds bounds])
    colormap bone
    
    subplot(9,5,26:35)
    imshow(cor_img,[-bounds bounds])
    colormap bone
    
    subplot(9,5,36:45)
    imshow(hor_img,[-bounds bounds])
    colormap bone
    
    print([savedir, label], '-dpng');
end

figure; 

BOLD_var = 0;
REJ_var = 0;
MID_var = 0;
IGN_var = 0;

numBOLD = size(accps,2);
numREJs = size(rejs,2);
numMIDS = size(mids,2);
numIGNS = size(igns,2);

for i = 1:numBOLD
    BOLD_var = BOLD_var+ imported_ctab(accps(i)+1,5);
end

for i = 1:numREJs
    REJ_var = REJ_var+ imported_ctab(rejs(i)+1,5);
end

for i = 1:numMIDS
    MID_var = MID_var+ imported_ctab((mids(i)+1),5);
end

for i = 1:numIGNS
    IGN_var = IGN_var + imported_ctab(igns(i)+1,5);
end

y = [BOLD_var, REJ_var, MID_var, IGN_var];

bar(1, y(1),  'facecolor', [0 .5 0]);
hold on
bar(2, y(2),  'facecolor', 'r');
bar(3, y(3),  'facecolor', 'm');
bar(4, y(4),  'facecolor', 'k');

boldLabel = strcat(num2str(numBOLD), ', BOLD');
noboldLabel = strcat(num2str(numREJs), ', non-BOLD');
r2Label = strcat(num2str(numMIDS), ', R2* Weighted');
ignLabel = strcat(num2str(numIGNS), ', Ignored');

labels = {boldLabel,noboldLabel,r2Label,ignLabel};
set(gca, 'XTick', 1:4, 'XTickLabel', labels);

%total_var = sum(y);
title(['% Exp. Var. of total ', total_var, ' PCA Variance ']);
ylabel('Variance Explained, %');

print([savedir, 'Var_exp'], '-dpng');

figure; 
plot(imported_ctab(:,2))
hold on
plot(imported_ctab(:,3))
legend('kappa', 'rho');
title(['Elbow, Kappa sorted vs Rho']);
ylabel('Value');
xlabel('Component Number');

print([savedir, 'Elbow_Graph_KappaVsRho'], '-dpng');

k_r = horzcat(imported_ctab(:,2),imported_ctab(:,3));
figure; 
scatter(k_r(:,1),k_r(:,2), imported_ctab(:,4)*100, color_table);
title(['Kappa vs Rho']);
ylabel('Rho');
xlabel('Kappa');

print([savedir, 'KappaVsRho'], '-dpng');

%%Lets make some tSNR figures as well, because why not. 
% At the moment there will be no filtering on these (highpass, etc)
% This could have ramifications for intepreting the data, but at the
% moment, this seems reasonable. 

cd(ted_dir);

tsoc_data = load_nii('ts_OC.nii');
tsoc_data = tsoc_data.img;
mean_tsoc = squeeze(mean(tsoc_data,4));
std_tsoc = squeeze(std(tsoc_data,0,4));
tsnr_tsoc = mean_tsoc./std_tsoc;

sag_img = [];
    for j = 1:sag_cuts:sag_slices
    sag_img = horzcat(sag_img, rot90(squeeze(tsnr_tsoc(j,:,:))));
    end
    
    cor_img = [];
    for j = 1:cor_cuts:cor_slices
    cor_img = horzcat(cor_img, rot90(squeeze(tsnr_tsoc(:,j,:))));
    end
    
    hor_img = [];
    for j = 1:hor_cuts:hor_slices
    hor_img = horzcat(hor_img, rot90(squeeze(tsnr_tsoc(:,:,j))));
    end

tsnr_range = max(max(max(tsnr_tsoc)))*.5;

figure; 
subplot(6,5,1:10)
imshow(sag_img,[0 tsnr_range])
colormap parula

subplot(6,5,11:20)
imshow(cor_img,[0 tsnr_range])
colormap parula

subplot(6,5,21:30)
imshow(hor_img,[0 tsnr_range])
colormap parula

h = colorbar; 
set(h, 'Position', [.07 .1 .03 .8150])

print([savedir, 'tsoc_tsnr'], '-dpng');

%Now for the Denoised Timeseries
medn_data = load_nii('dn_ts_OC.nii');
medn_data = medn_data.img;
mean_medn = squeeze(mean(medn_data,4));
std_medn = squeeze(std(medn_data,0,4));
tsnr_medn = mean_medn./std_medn;

sag_img = [];
    for j = 1:sag_cuts:sag_slices
    sag_img = horzcat(sag_img, rot90(squeeze(tsnr_medn(j,:,:))));
    end
    
    cor_img = [];
    for j = 1:cor_cuts:cor_slices
    cor_img = horzcat(cor_img, rot90(squeeze(tsnr_medn(:,j,:))));
    end
    
    hor_img = [];
    for j = 1:hor_cuts:hor_slices
    hor_img = horzcat(hor_img, rot90(squeeze(tsnr_medn(:,:,j))));
    end

tsnr_range = max(max(max(tsnr_medn)))*.5;

figure; 
subplot(6,5,1:10)
imshow(sag_img,[0 tsnr_range])
colormap parula

subplot(6,5,11:20)
imshow(cor_img,[0 tsnr_range])
colormap parula

subplot(6,5,21:30)
imshow(hor_img,[0 tsnr_range])
colormap parula

h = colorbar; 
set(h, 'Position', [.07 .1 .03 .8150])

print([savedir, 'medn_tsnr'], '-dpng');

%and one last time for the ratio between the two. 


tsnr_ratio = tsnr_medn./tsnr_tsoc;
sag_img = [];
    for j = 1:sag_cuts:sag_slices
    sag_img = horzcat(sag_img, rot90(squeeze(tsnr_ratio(j,:,:))));
    end
    
    cor_img = [];
    for j = 1:cor_cuts:cor_slices
    cor_img = horzcat(cor_img, rot90(squeeze(tsnr_ratio(:,j,:))));
    end
    
    hor_img = [];
    for j = 1:hor_cuts:hor_slices
    hor_img = horzcat(hor_img, rot90(squeeze(tsnr_ratio(:,:,j))));
    end

tsnr_range = max(max(max(tsnr_ratio)))*.3;

figure; 
subplot(6,5,1:10)
imshow(sag_img,[1 tsnr_range])
colormap parula

subplot(6,5,11:20)
imshow(cor_img,[1 tsnr_range])
colormap parula

subplot(6,5,21:30)
imshow(hor_img,[1 tsnr_range])
colormap parula

h = colorbar; 
set(h, 'Position', [.07 .1 .03 .8150])

print([savedir, 'tsnr_ratio'], '-dpng');

