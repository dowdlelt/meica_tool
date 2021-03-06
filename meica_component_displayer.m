function meica_component_displayer(tr, savedir)

% tr is the repititon time of the data.
% savedir is the TED output folder that contains everything tedana produced
if exist('savedir','var')

else
    savedir = spm_select(1,'dir','Select the MEICA output folder...');
end
%Make sure to call it with a tr now.
%and the figures will be saved in new folder here.

cd(savedir);
try copyfile ../../dfile.r01.1D ./;
catch
    fprintf('dfile not found')

end

% Custom color map
[path,name,ext] = fileparts(which('meica_component_displayer'));
custom_cmap_loc = fullfile(path,'custom_map.mat');
load(custom_cmap_loc);


%Prepare motion calc and loading.
cfg.motionparam = 'dfile.r01.1D'; %output from MEICA, organized as: roll pitch yaw dS  dL  dP
cfg.prepro_suite = 'meica';
cfg.radius = 50;

Motion = 1;

try raw_motion = load(cfg.motionparam); %Load up motion for plotting.
catch no_motion
    fprintf('no motion file found, ignoring motionrelated plots\n');
    Motion = 0;
end


% ted_dir = strcat(savedir, '/TED.r01');
ted_dir = savedir;
cd(ted_dir);

%Get the ICA component timecourses.
timecourses_data = load('meica_mix.1D');
ctab = 'comp_table.txt';


if Motion
%Create a combination of motion estimates and ICA comps.
total_timecourses = horzcat(raw_motion, timecourses_data);
%Use this to create a correlation matrix
corr_mat = corr(total_timecourses);
end

fid = fopen(ctab); %Open the comp table that is within the TED folder.
tline = fgetl(fid);
imported_ctab = [];

while ischar(tline)
    %disp(tline)
    num_check = str2num(tline(1)); %Is the first character a number

    if isempty(num_check) %Want to make sure it is not []
        num_check = 'not a number'; %not a number ignore it
    end

    if strfind(tline, '#ACC')
        accp_list = tline;
    elseif strfind(tline, '#REJ')
        rej_list = tline;
    elseif strfind(tline, '#MID')
        mid_list = tline;
    elseif strfind(tline, '#IGN')
        ign_list = tline;
    elseif strfind(tline, '(VEx)')
        total_var = tline;
    elseif isnumeric(num_check) %if first character is number - its a component
        imported_ctab = vertcat(imported_ctab, str2num(tline));%Build the component table
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

%For the imported_ctab, the columns are
%Comp#  Kappa   Rho   Variance   Normed_variance

cd(ted_dir); %Go into the output directory

try
    all_betas = load_nii('betas_OC.nii.gz');
catch
    warning('uncompressed niis')
    all_betas = load_nii('betas_OC.nii');
end

all_betas = all_betas.img;

cd(savedir); %Go back out into the main directory.
if Motion
%%Making the motion plots, including framewise displacement

figure('visible','off', 'windowstyle', 'normal');
grid on; grid minor;
[fwd,~]=bramila_framewiseDisplacement(cfg); %calculate FD using script

x_axis = size(fwd,1); %Get the number of timepoints

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
end
mkdir('component_plots');
savedir = [savedir, '/component_plots/'];

%plots are created here.
[~, titl, ~] = fileparts(cfg.motionparam);
print([savedir, titl], '-dpng');

x_axis = size(timecourses_data,1);
comp_number = size(timecourses_data,2);

num_cuts = 10; %This is used to creating multi-planar images

color_table = zeros(comp_number,3);

fprintf('Making component plots.\n');
count = 0; %for displaying things every so often

figure('visible','off', 'windowstyle','normal');

n_colors = 256;
colormaptype = 'parula';
cmap_fh = str2func(colormaptype);

for i = 1:comp_number
%    cmap = cmap_fh(n_colors);
%     for v = 1:(n_colors/2)
%         cmap(v,:) = [1 1 1];
%     end
    cmap = custom_map;

    if mod(count, 10) == 0
        fprintf('.'); %print progress every 10 components, less annoying.
    end
    count = count +1;

    upper = max(max(max(all_betas(:,:,:,i))));
    lower = min(min(min(all_betas(:,:,:,i))));
    bounds = (0.3 * max([abs(upper), abs(lower)]));

    kappa = imported_ctab(i,2);
    rho = imported_ctab(i,3);
    variance_explained = imported_ctab(i,4);

    %% Prints the graphs so that the outcome can be seen
    subplot(12,5,1:15)
    length_of_time = size(timecourses_data(:,i),1);

    if any(accps == (i-1))
        plot(timecourses_data(:,i), 'Color',[0 .5 0]); %Green for BOLD like
        color_table(i,:) = [0 1 0];
    elseif any(rejs == (i-1))
        plot(timecourses_data(:,i), 'r'); %Red for Rejected non BOLD
        color_table(i,:) = [1 0 0];
    elseif any(mids == (i-1))
        plot(timecourses_data(:,i), 'm'); %Magenta for R2* weighted artifacts
        color_table(i,:) = [1 0 1];
    elseif any(igns == (i-1))
        plot(timecourses_data(:,i), 'k'); %Black for Ignored components
        color_table(i,:) = [0 0 0];
    end
    %%

    axis([0 x_axis min(timecourses_data(:,i)) max(timecourses_data(:,i))]);
    title(strcat('Component:', num2str(i), ', on ctab: ', num2str(i-1), ', kappa: ', num2str(kappa,3), ', rho: ', num2str(rho,3), ', variance: ', num2str(variance_explained,4))); grid on;
    label = strcat('Component_', num2str(i), '_on_ctab_', num2str(i-1));

    current_image = squeeze(all_betas(:,:,:,i));

    [sag_img, cor_img, hor_img] = three_cut_maker(current_image,num_cuts);

    subplot(12,5,16:25)
    imshow(sag_img,[-bounds bounds],'Colormap',cmap)

    subplot(12,5,26:35)
    imshow(cor_img,[-bounds bounds],'Colormap',cmap)

    subplot(12,5,36:45)
    imshow(hor_img,[-bounds bounds],'Colormap',cmap)

    %Now adding in the fourier transform.
    Fs = 1/tr;            % Sampling frequency
    T = 1/Fs;             % Sampling period
    L = (length_of_time*T);             % Length of signal

    Y = fft(timecourses_data(:,i));
    P2 = abs(Y/L);
    P1 = P2(1:floor(L/2+1));
    P1(2:end-1) = 2*P1(2:end-1);

    f = Fs*(0:(L/2))/L;
    subplot(12,5,52:60)
    plot(f,P1)
    title('Single-Sided Amplitude Spectrum of X(t)')
    xlabel('f (Hz)')
    ylabel('|P1(f)|')

    print([savedir, label], '-dpng');
end
%%
fprintf('\nCollecting Explained Variance...');

figure('visible','off', 'windowstyle','normal');

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

title(['% Exp. Var. of total ', total_var, ' PCA Variance ']);
ylabel('Variance Explained, %');

print([savedir, 'Var_exp'], '-dpng');

%%
fprintf('\nShowing Elbow of Kappa, with Rho...');

figure('visible','off', 'windowstyle','normal');
plot(imported_ctab(:,2))
hold on
plot(imported_ctab(:,3))
legend('kappa', 'rho');
title('Elbow, Kappa sorted vs Rho');
ylabel('Value');
xlabel('Component Number');

print([savedir, 'Elbow_Graph_KappaVsRho'], '-dpng');
%%
fprintf('\nScattering Kappa vs Rho ');

k_r = horzcat(imported_ctab(:,2),imported_ctab(:,3));
figure('visible','off', 'windowstyle','normal');
scatter(k_r(:,1),k_r(:,2), imported_ctab(:,4)*100, color_table);
title('Kappa vs Rho');
ylabel('Rho');
xlabel('Kappa');

print([savedir, 'KappaVsRho'], '-dpng');

%%
%Lets make some tSNR figures as well, because why not.
% At the moment there will be no filtering on these (highpass, etc)
% This could have ramifications for intepreting the data, but at the
% moment, this seems reasonable.

cd(ted_dir);

try
    base_img = load_nii('t2sv.nii.gz');
catch
    base_img = load_nii('t2sv.nii');
end
%This is a one frame nifti that we can use to make nifti versions
% of all the TSNR figures.

fprintf('\nCalculating TSNR figures...');

%%
%Calculated TSNR denoised Timeseries
try
    tsnr_medn = tsnr_creator('dn_ts_OC.nii.gz');
catch
    tsnr_medn = tsnr_creator('dn_ts_OC.nii');
end

[sag_img, cor_img, hor_img] = three_cut_maker(tsnr_medn,num_cuts);

tsnr_max_MEDN = max(max(max(tsnr_medn)));

tsnr_range = tsnr_max_MEDN*0.8;

figure('visible','off', 'windowstyle','normal');
subplot(8,5,1:10)
imshow(sag_img,[0 tsnr_range])
title('TSNR of MEDN timeseries');
colormap parula

subplot(8,5,11:20)
imshow(cor_img,[0 tsnr_range])
colormap parula

subplot(8,5,21:30)
imshow(hor_img,[0 tsnr_range])
colormap parula

h = colorbar;
set(h, 'Position', [.08 .35 .03 .55])

subplot(8,5,31:40)
histogram(reshape(tsnr_medn, [],1),100);
xlim([0 tsnr_max_MEDN]);

print([savedir, 'tsnr_medn'], '-dpng');

base_img.img = tsnr_medn;
save_nii(base_img, [savedir, 'medn_tsnr.nii']);

%%
%Calculate TSNR of TSOC
if size(dir('tsoc_orig.nii*'),1) >= 1 %older tedana versions
    try
        tsnr_tsoc = tsnr_creator('tsoc_orig.nii.gz');
    catch
        tsnr_tsoc = tsnr_creator('tsoc_orig.nii');
    end
elseif size(dir('ts_OC.nii*'),1) >= 1 %newer tedana versions
    try
        tsnr_tsoc = tsnr_creator('ts_OC.nii.gz');
    catch
        tsnr_tsoc = tsnr_creator('ts_OC.nii');
    end
else
    disp('Optimally combined file not found');
end

[sag_img, cor_img, hor_img] = three_cut_maker(tsnr_tsoc,num_cuts);

figure('visible','off', 'windowstyle','normal');
subplot(8,5,1:10)
imshow(sag_img,[0 tsnr_range])

title('TSNR of TSOC timeseries');
colormap parula

subplot(8,5,11:20)
imshow(cor_img,[0 tsnr_range])
colormap parula

subplot(8,5,21:30)
imshow(hor_img,[0 tsnr_range])
colormap parula

h = colorbar;
set(h, 'Position', [.08 .35 .03 .55])

subplot(8,5,31:40)
histogram(reshape(tsnr_tsoc, [],1),100);
xlim([0 tsnr_max_MEDN]);

print([savedir, 'tsnr_oc'], '-dpng');

base_img.img = tsnr_tsoc;
save_nii(base_img, [savedir, 'tsoc_tsnr.nii']);
%%
%Calculated TSNR of second echo
%This is typically around 30
%So its tSNR is close to what a conventional aquisition would be.
cd ..

tsnr_e2 = tsnr_creator('e2_in.nii.gz');

[sag_img, cor_img, hor_img] = three_cut_maker(tsnr_e2,num_cuts);

figure('visible','off', 'windowstyle','normal');
subplot(8,5,1:10)
imshow(sag_img,[0 tsnr_range])
title('TSNR of Second Echo');
colormap parula

subplot(8,5,11:20)
imshow(cor_img,[0 tsnr_range])
colormap parula

subplot(8,5,21:30)
imshow(hor_img,[0 tsnr_range])
colormap parula

h = colorbar;
set(h, 'Position', [.08 .35 .03 .55])

subplot(8,5,31:40)
histogram(reshape(tsnr_e2, [],1),100);
xlim([0 tsnr_max_MEDN]);

print([savedir, 'tsnr_2nd_echo'], '-dpng');

base_img.img = tsnr_e2;
save_nii(base_img, [savedir, 'e2_tsnr.nii']);
%%
%Ratio of Denoised vs TSOC TSNR

fprintf('\nCalculating TSNR ratio');
tsnr_ratio = tsnr_medn./tsnr_tsoc;

[sag_img, cor_img, hor_img] = three_cut_maker(tsnr_ratio,num_cuts);

figure('visible','off', 'windowstyle','normal');
subplot(8,5,1:10)
imshow(sag_img,[0 5])
title('TSNR Ratio, MEDN vs TSOC timeseries');
colormap parula

subplot(8,5,11:20)
imshow(cor_img,[0 5])
colormap parula

subplot(8,5,21:30)
imshow(hor_img,[0 5])
colormap parula

h = colorbar;
set(h, 'Position', [.08 .35 .03 .55])

subplot(8,5,31:40)
c = histogram(reshape(tsnr_ratio, [],1),25);

c.BinLimits = [0 6];
c.NumBins = 30;

print([savedir, 'tsnr_ratio_medn_tsoc'], '-dpng');

base_img.img = tsnr_ratio;
save_nii(base_img, [savedir, 'tsnr_ratio_medn_tsoc.nii']);

cd ..
%%
%TSNR ratio from optimal combo vs single echo

fprintf('\nCalculating TSNR ratio for TSOC vs E2');
tsnr_ratio = tsnr_tsoc./tsnr_e2;

[sag_img, cor_img, hor_img] = three_cut_maker(tsnr_ratio,num_cuts);

figure('visible','off', 'windowstyle','normal');
subplot(8,5,1:10)
imshow(sag_img,[0 5])
title('TSNR of TSOC vs 2nd Echo');
colormap parula

subplot(8,5,11:20)
imshow(cor_img,[0 5])
colormap parula

subplot(8,5,21:30)
imshow(hor_img,[0 5])
colormap parula

h = colorbar;
set(h, 'Position', [.08 .35 .03 .55])

subplot(8,5,31:40)
c= histogram(reshape(tsnr_ratio, [],1),25);
c.BinLimits = [0 6];
c.NumBins = 30;

print([savedir, 'tsnr_ratio_oc_vs_e2'], '-dpng');
cd ..

%%
%TSNR ratio for denoising vs single echo
fprintf('\nCalculating TSNR ratio for MEDN vs E2');

tsnr_ratio = tsnr_medn./tsnr_e2;

[sag_img, cor_img, hor_img] = three_cut_maker(tsnr_ratio,num_cuts);

figure('visible','off', 'windowstyle','normal');
subplot(8,5,1:10)
imshow(sag_img,[0 5])
title('TSNR of MEDN vs 2nd Echo');
colormap parula

subplot(8,5,11:20)
imshow(cor_img,[0 5])
colormap parula

subplot(8,5,21:30)
imshow(hor_img,[0 5])
colormap parula

h = colorbar;
set(h, 'Position', [.08 .35 .03 .55])

subplot(8,5,31:40)
c = histogram(reshape(tsnr_ratio, [],1),25);
c.BinLimits = [0 6];
c.NumBins = 30;

print([savedir, 'tsnr_ratio_medn_vs_e2'], '-dpng');

base_img.img = tsnr_ratio;
save_nii(base_img, [savedir, 'tsnr_ratio_medn_vs_e2.nii']);

%Creating a covariance matrix
fprintf('\nPlotting correlations');
figure('visible','off', 'windowstyle','normal');
imagesc(corr_mat, [-1 1]);
title('Correlation of Motion and ICA comps.');
colorbar
colormap parula
print([savedir, 'Covariance_matrix'], '-dpng');

%%
%Creating R^2 matrix
fprintf('\nPlotting correlations...squared\n');
figure('visible','off', 'windowstyle','normal');
imagesc((corr_mat.^2), [0 1]);
title('Correlation Coefficient of Motion and ICA comps.');
colorbar
colormap parula
print([savedir, 'R_Squared_matrix'], '-dpng');
