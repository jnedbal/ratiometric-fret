function FRETeff
global files
global channels
global fs

% six samples:
%   1. Background         - cells transfected with empty plasmid pCDNA MH
%   2. Donor only         - cells transfected with GFP
%   3. Acceptor only      - cells transfected with RFP
%   4. Donor and Acceptor - cells transfected with 7aa FRET rulers
%   5. Donor and Acceptor - cells transfected with 19aa FRET rulers
%   6. Donor and Acceptor - cells transfected with 32aa FRET rulers

% Three detection channels
%   I1. donor emission (488 nm ex., 510 ± 7.5 nm em.)
%   I2. sensitized acceptor emission (488 nm ex., 585 ± 21 nm em.)
%   I1. acceptor emission (561 nm ex., 582 ± 7.5 nm em.)

%% Set global font size
fs = 20;

%% Set files
% prefix of the file common to all
prefix = 'Voltage-2_';
% file extension
extension = '.fcs';
% file name identifying parts
names = {'PC-DNA', 'gfp', 'rfp', '7aa', '19aa', '32aa'};
% type of data, match with the names above to give them meaning:
%   0:  background
%   1:  acceptor
%   2:  donor
%   3:  donor-acceptor
types = [0 2 1 3 3 3];
% names of flow cytometer channels to use
channels.name.ssc = 'SSC-A';         % side scatter amplitude channel
channels.name.sscw = 'SSC-W';         % side scatter width channel
channels.name.fsc = 'FSC-A';         % forward scatter channel
channels.name.donor = 'FITC-A';      % donor fluorescence
channels.name.acceptor = 'PE-A';     % acceptor fluorescence
channels.name.transfer = 'PerCP-A';  % sensitized acceptor fluorescence

% For log graphs, add an offset to make sure that we don't have negative
% values
logOff = 80;
% Transfection gates value
transGate = [180, 180];

%% Create folder for images
if ~(exist([pwd filesep 'images'], 'dir') == 7)
    mkdir('images');
end

%% Load data
for i = 1 : numel(names)
    files(i).name = strcat(prefix, names{i}, extension);
    [files(i).fcsdat, ...
        files(i).fcshdr, ...
        files(i).fcsdatscaled, ...
        files(i).fcsdat_comp] = ...
                                fca_readfcs(files(i).name);
end

%% Get indices of the channels that are needed
% assume that all samples have the same number and order of channels.
in = types == 0;    % background channel
chnames = {files(in).fcshdr.par.name};
chdesc = fields(channels.name);
for i = 1 : numel(chdesc)
    channels.index.(chdesc{i}) = ...
        cellfun(@(x) strcmpi(x, channels.name.(chdesc{i})), chnames);
end





% *********************************************************************** %
% *                    Forward - Side Scatter Plots                     * %
% *********************************************************************** %



%% Draw the density dot plots for forward and side scatter 
% create a forward-side scatter gate. It contains a series of X, Y
% coordinates of points that define the gate polygon
par.gate = [ 60000,  45000; ...
            100000,  45000; ...
            140000,  75000; ...
            180000, 133000; ...
            147000, 165000; ...
             92000, 137000; ...
             60000,  86000];
par.fiden = 'SSC-FSC';      % filename identifier of the output image
par.xaxis = 'Forward Scatter (FSC)';    % X-axis label
par.yaxis = 'Side Scatter (SSC)';       % Y-axis label
par.title = 'Forward-Side';             % graph title
par.colormap = 'autumn';                % colormap of the gated dots

for i = 1 : numel(names)
    % set parameters for forward-side scatter plots
    par.type = types(i);
    if par.type == 3
        par.name = names{i};
    else
        par.name = '';
    end
    % filename index
    par.findex = i;
    % Forward scatter
    par.xdata = files(i).fcsdat(:, channels.index.fsc);
    % Side scatter
    par.ydata = files(i).fcsdat(:, channels.index.ssc);
    % axis limits
    par.xlim = [0, files(i).fcshdr.par(channels.index.fsc).range];
    par.ylim = [0, files(i).fcshdr.par(channels.index.ssc).range];
    % filename to print on the graph
    par.fname = files(i).fcshdr.filename;

    % create forward/side scatter plots
    gate = scatter_plot(par);

    % FSC-SSC gate for the rest of the experiments
    files(i).FSgate = gate;
end
clear par    % clear to prevent carry-over to next graphs





% *********************************************************************** %
% *                          Singlet Cell Gate                          * %
% *********************************************************************** %



u = i + 1;          % graph index
%% Draw the density dot plots for forward and side scatter 
% create a forward-side scatter gate. It contains a series of X, Y
% coordinates of points that define the gate polygon
par.gate = [ 70000,  45000; ...
            100000,  45000; ...
            120000, 165000; ...
             90000, 165000];
par.fiden = 'SSCW-SSCA';      % filename identifier of the output image
par.xaxis = 'Side Scatter Width (SSC-W)';        % X-axis label
par.yaxis = 'Side Scatter Amplitude (SSC-A)';    % Y-axis label
par.title = 'Side Width-Side Amplitude';         % graph title
par.colormap = 'autumn';                % colormap of the gated dots

for i = 1 : numel(names)
    % set parameters for forward-side scatter plots
    par.type = types(i);
    if par.type == 3
        par.name = names{i};
    else
        par.name = '';
    end
    % filename index
    par.findex = u; u = u + 1;
    % Forward scatter
    par.xdata = files(i).fcsdat(:, channels.index.sscw);
    % Side scatter
    par.ydata = files(i).fcsdat(:, channels.index.ssc);
    % axis limits
    par.xlim = [0, files(i).fcshdr.par(channels.index.sscw).range];
    par.ylim = [0, files(i).fcshdr.par(channels.index.ssc).range];
    % filename to print on the graph
    par.fname = files(i).fcshdr.filename;

    % create forward/side scatter plots
    gate = scatter_plot(par);

    % SSC-W SSC-A gate for the rest of the experiments
    files(i).SWSAgate = gate;

    % Combined scatter gate
    files(i).SCATgate = gate & files(i).FSgate;
end
clear par    % clear to prevent carry-over to next graphs





% *********************************************************************** %
% *                       Background Subtraction                        * %
% *********************************************************************** %



%% Calculate the median of the three channels of the background
% get the background channel
type = types == 0;  % background channel
channels.Bg.donor = ...
    median(files(type).fcsdat(files(type).SCATgate, channels.index.donor));
channels.Bg.transfer = ...
    median(files(type).fcsdat(files(type).SCATgate,channels.index.transfer));
channels.Bg.acceptor = ...
    median(files(type).fcsdat(files(type).SCATgate,channels.index.acceptor));





% *********************************************************************** %
% *                       Gate Transfected Cells                        * %
% *********************************************************************** %



tC = [0 0 0 0];     % type counter
axLab = {'Donor', 'Transfer', 'Acceptor'};  % Graph and axis title

for i = 1 : numel(types)
    %% get the type of the channel
    switch types(i)
        case 0
            type = 'Bg';
        case 1
            type = 'acceptor';
        case 2
            type = 'donor';
        case 3
            type = 'transfer';
    end

    %% increment the type counter
    tC(types(i) + 1) = tC(types(i) + 1) + 1;
    j = tC(types(i) + 1);

    %% Calculate the median of the three channels of the 
    channels.(type)(j).F1 = ...
        files(i).fcsdat(files(i).SCATgate, channels.index.donor) - ...
        channels.Bg(1).donor + logOff;
    channels.(type)(j).F2 = ...
        files(i).fcsdat(files(i).SCATgate, channels.index.transfer) - ...
        channels.Bg(1).transfer + logOff;
    channels.(type)(j).F3 = ...
        files(i).fcsdat(files(i).SCATgate, channels.index.acceptor) - ...
        channels.Bg(1).acceptor + logOff;

    %% Create a positive gate
    channels.(type)(j).posGate = true(size(channels.(type)(j).F1));

    for Fs = [1, 1; 2, 3];
        %% Plot the background corrected intensities
        par.quadrant = transGate;
        par.mergeq = [1, 2, 2, 2];
        par.dscatter = {'logxlogy', true};
        % filename identifier of the output image
        par.fiden = sprintf('F%d - F%d', Fs);
        % X-axis label
        xLab = axLab{Fs(1)};
        par.xaxis = sprintf('%s Channel Intensity (F%d)', xLab, Fs(1));
        % Y-axis label
        yLab = axLab{Fs(2)};
        par.yaxis = sprintf('%s Channel Intensity (F%d)', yLab, Fs(2));
        % graph title
        par.title = sprintf('%s-%s Intensity', yLab, xLab);

        % set parameters for the scatter plot
        % background (0), acceptor (1), donor (2), transfer(3)
        par.type = types(i);
        if par.type == 3
            par.name = names{i};
        else
            par.name = '';
        end

        % filename index
        par.findex = u; u = u + 1;
        % X Channel Intensity
        par.xdata = channels.(type)(j).(sprintf('F%d', Fs(1)));
        % Y Channel Intensity
        par.ydata = channels.(type)(j).(sprintf('F%d', Fs(2)));
        % axis limits
        par.xlim = ...
            [10, files(i).fcshdr.par(channels.index.(lower(xLab))).range];
        par.ylim = ...
            [10, files(i).fcshdr.par(channels.index.(lower(yLab))).range];
        % colormap
        par.colormap = 'autumn';
        % filename to print on the graph
        par.fname = files(i).fcshdr.filename;

        % create the scatter plot
        gate = scatter_plot(par);

        % positive gate
        channels.(type)(j).(sprintf('posF%dF%dgate', Fs)) = gate;
        channels.(type)(j).posGate = channels.(type)(j).posGate & ...
                                        ~gate(:, 1);
        clear par    % clear to prevent carry-over to next graphs
    end
end





% *********************************************************************** %
% *                     Process the gated data                          * %
% *********************************************************************** %



tC = [0 0 0 0];     % type counter
for i = 1 : numel(types)
    %% No need to plot background sample data
    if types(i) == 0
        continue
    end
    %% get the type of the channel
    switch types(i)
        case 0
            type = 'Bg';
        case 1
            type = 'acceptor';
        case 2
            type = 'donor';
        case 3
            type = 'transfer';
    end

    %% increment the type counter
    tC(types(i) + 1) = tC(types(i) + 1) + 1;
    j = tC(types(i) + 1);

    %% Create a non-negative gate
    channels.(type)(j).nonnegGate = ...
        true(sum(channels.(type)(j).posGate), 1);

    for Fs = [1, 1, 2; 2, 3, 3];
        if Fs(1) == 2 && types(i) ~= 3
            continue
        end
        %% Plot the background corrected Background intensities F1 and F2
        par.quadrant = [logOff logOff];
        par.mergeq = [1, 2, 2, 4];
        par.dscatter = {'logxlogy', true};
        % filename identifier of the output image
        par.fiden = sprintf('F%d - F%d', Fs);
        % X-axis label
        xLab = axLab{Fs(1)};
        par.xaxis = sprintf('%s Channel Intensity (F%d)', xLab, Fs(1));
        % Y-axis label
        yLab = axLab{Fs(2)};
        par.yaxis = sprintf('%s Channel Intensity (F%d)', yLab, Fs(2));
        % graph title
        par.title = sprintf('%s-%s Intensity', yLab, xLab);

        % set parameters for scatter plots
        % background (0), acceptor (1), donor (2), transfer(3)
        par.type = types(i);
        if par.type == 3
            par.name = names{i};
        else
            par.name = '';
        end

        % filename index
        par.findex = u; u = u + 1;
        % X Channel Intensity
        par.xdata = channels.(type)(j).(sprintf('F%d', Fs(1)))(channels.(type)(j).posGate);
        % Y Channel Intensity
        par.ydata = channels.(type)(j).(sprintf('F%d', Fs(2)))(channels.(type)(j).posGate);
        % axis limits
        par.xlim = ...
            [10, files(i).fcshdr.par(channels.index.(lower(xLab))).range];
        par.ylim = ...
            [10, files(i).fcshdr.par(channels.index.(lower(yLab))).range];
        % colormap
        par.colormap = copper(128);
        par.colormap = vertcat({'autumn'}, ...
                               repmat({par.colormap(65 : end, :)}, 2, 1));
        % filename to print on the graph
        par.fname = files(i).fcshdr.filename;

        % create scatter plot
        gate = scatter_plot(par);

        % non-negative gate
        channels.(type)(j).(sprintf('nonnegF%dF%dgate', Fs)) = gate;
        channels.(type)(j).nonnegGate = channels.(type)(j).nonnegGate & ...
                                        gate(:, 4);

        clear par    % clear to prevent carry-over to next graphs

    end
end





% *********************************************************************** %
% *                   Acceptor Channel Analysis                         * %
% *********************************************************************** %



%% Calculate the S2 and S4 bleed-through values from the ACCEPTOR sample

% S2: bleed-though of the direct acceptor emission into
%     sensitized (transfer) channel.
channels.S2 = [];
% S4: bleed-though of the direct acceptor emission into
%     donor channel.
channels.S4 = [];
% MA: median of the sensitized ACCEPTOR emission
channels.MA = [];
for i = 1 : numel(channels.acceptor)
    F1 = channels.acceptor(i).F1 - logOff;
    F1 = F1(channels.acceptor(i).posGate);
    F1 = F1(channels.acceptor(i).nonnegGate);
    F2 = channels.acceptor(i).F2 - logOff;
    F2 = F2(channels.acceptor(i).posGate);
    F2 = F2(channels.acceptor(i).nonnegGate);
    F3 = channels.acceptor(i).F3 - logOff;
    F3 = F3(channels.acceptor(i).posGate);
    F3 = F3(channels.acceptor(i).nonnegGate);

    channels.S2 = vertcat(channels.S2, F2 ./ F3);
    channels.S4 = vertcat(channels.S4, F1 ./ F3);
    channels.MA = vertcat(channels.MA, F2);
end
channels.MA = median(channels.MA);





% *********************************************************************** %
% *                     Donor Channel Analysis                          * %
% *********************************************************************** %



%% Calculate the S1 and S3 bleed-through values from the DONOR sample

% S1: bleed-though of the donor emission into
%     sensitized (transfer) channel.
channels.S1 = [];
% S3: bleed-though of the donor emission into the
%     direct excitation acceptor channel.
channels.S3 = [];
% MD: median of the DONOR emission
channels.MD = [];
for i = 1 : numel(channels.donor)
    F1 = channels.donor(i).F1 - logOff;
    F1 = F1(channels.donor(i).posGate);
    F1 = F1(channels.donor(i).nonnegGate);
    F2 = channels.donor(i).F2 - logOff;
    F2 = F2(channels.donor(i).posGate);
    F2 = F2(channels.donor(i).nonnegGate);
    F3 = channels.donor(i).F3 - logOff;
    F3 = F3(channels.donor(i).posGate);
    F3 = F3(channels.donor(i).nonnegGate);

    channels.S1 = vertcat(channels.S1, F2 ./ F1);
    channels.S3 = vertcat(channels.S3, F3 ./ F1);
    channels.MD = vertcat(channels.MD, F1);
end
channels.MD = median(channels.MD);





% *********************************************************************** %
% *                        Alpha Coefficient                            * %
% *********************************************************************** %



%% Calculate the alpha value
% alpha = (MA / MD) * (LD / LA) * (ED / EA)
%
%   alpha is the coefficient
%   MA is the sensitized acceptor emission
%   MD is the donor emission
%   LD is labeling ratio of the donor
%   LA is labeling ratio of the acceptor
%   ED is the donor molar excitation coefficient at donor wavelength
%   EA is the acceptor molar excitation coefficient at donor wavelength
channels.LD = 1;
channels.LA = 1;
%http://www.tsienlab.ucsd.edu/Documents/REF%20-%20Fluorophore%20Spectra.xls
channels.ED = 0.998 * 55000;
% Biophys J Marion Peter
channels.EA = 13048;
channels.EA = 0.12 * 44000;

% Calculate alpha
channels.alpha = channels.MA / channels.MD * ...
                 channels.LD / channels.LA * ...
                 channels.ED / channels.EA;





% *********************************************************************** %
% *                    S Bleed-Through Correction                       * %
% *********************************************************************** %



%% Plot histograms of the S bleed-through corrections
for i = 1 : 4
    par.gate = [0, 0.3];            % gate for S
    % filename identifier of the output image
    par.fiden = sprintf('S%d', i);
    par.yaxis = 'Frequency';        % Y-axis label
    switch i
        case 1
            % X-axis label
            par.xaxis = 'S1 = F2 / F1';
            % graph title
            par.title = 'Donor to Sensitized Emission Bleed-Through';
            % donor
            par.type = 2;
        case 2
            % X-axis label
            par.xaxis = 'S2 = F2 / F3';
            % graph title
            par.title = 'Acceptor to Sensitized Emission Bleed-Through';
            % acceptor
            par.type = 1;
        case 3
            % X-axis label
            par.xaxis = 'S3 = F3 / F1';
            % graph title
            par.title = 'Donor to Acceptor Emission Bleed-Through';
            % donor
            par.type = 2;
        case 4
            % X-axis label
            par.xaxis = 'S4 = F1 / F3';
            % graph title
            par.title = 'Acceptor to Donor Emission Bleed-Through';
            % acceptor
            par.type = 1;
    end
    par.color = [0 0.75 0];         % Gate Color
    % set parameters for F3-F1 scatter plots

    par.name = '';

    % filename index
    par.findex = u; u = u + 1;


    % input data
    par.data = channels.(sprintf('S%d', i));
    % histogram bins (either a vector of bins, or number of bins)
    par.xbins = 0 : 0.01 : 1;
    par.xlim = [0 1];   % X-axis limits
    % filename to print on the graph
    par.fname = files(types == par.type).fcshdr.filename;

    % create histogram
    gate = histogram(par);

    % Calculate the S value median
    channels.(sprintf('S%dmed', i)) = median(par.data(gate));

    clear par    % clear to prevent carry-over to next graphs
end





% *********************************************************************** %
% *                      R_1 and R_F Calculation                        * %
% *********************************************************************** %



% R_F = F_AD/F_A where F_AD = F2-I1*S1 and F_A ~ F3*S2
% R_I = (F2-S2F3)/((1-S2S3/S1)F1)-S1

S1m = channels.S1med;
S2m = channels.S2med;
S3m = channels.S3med;

close all

X = zeros(1, numel(channels.transfer));
Y = X;

for i = 1 : numel(X)
    F1 = channels.transfer(i).F1 - logOff;
    F1 = F1(channels.transfer(i).posGate);
    F1 = F1(channels.transfer(i).nonnegGate);
    F2 = channels.transfer(i).F2 - logOff;
    F2 = F2(channels.transfer(i).posGate);
    F2 = F2(channels.transfer(i).nonnegGate);
    F3 = channels.transfer(i).F3 - logOff;
    F3 = F3(channels.transfer(i).posGate);
    F3 = F3(channels.transfer(i).nonnegGate);
    F_A = F3 * S2m;
    F_AD = F2 - F1 * S1m;
    channels.transfer(i).RF = F_AD ./ F_A;
    channels.transfer(i).RI = (F2 - S2m * F3) ./ ((1 - S2m * S3m / S1m) * F1) - S1m;
    X(i) = median(1 ./ channels.transfer(i).RI);
    Y(i) = median(1 ./ (channels.transfer(i).RF - 1));
end

% Fit a line through X/Y plot
f1 = figure('Units', 'Pixels', 'Position', [100 100 600 600]);
axes('Units', 'Pixels', 'Position', [110 70 450 450], 'FontSize', fs);
plot(X, Y, 'ko', 'MarkerSize', 10)
hold on
leg = {'7aa', '19aa', '32aa'};  % legend

for i = 1 : numel(X)
    text(X(i), Y(i), [leg{i} '  '], 'FontSize', fs, ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom')
end
p = polyfit(X, Y, 1);
plot([0 ceil(X(end))], [0 ceil(X(end))] * p(1) + p(2), 'k-', 'LineWidth', 2)
set(gca, 'YLim', [0, max(get(gca, 'YLim'))])
yt = get(gca, 'YTick');
xt = get(gca, 'XTick');
ytl = get(gca, 'YTickLabel');
ytl(abs(yt - p(2)) == min(abs(yt - p(2))), :) = repmat(' ', 1, size(ytl, 2));
set(gca, 'YTicklabel', ytl);
plot([0 ceil(X(end)) / 5], p(2) * [1 1], 'k:', 'LineWidth', 2)
text(xt(end) / 5, p(2), sprintf(' %4.3f', p(2)), 'FontSize', fs)
text(xt(end) / 5, xt(end) / 5 * p(1) + p(2), ...
     sprintf(' %4.3f', p(1)), 'FontSize', fs, 'Rotation', ...
     atan(p(1)) * 180 / pi * diff(get(gca, 'XLim')) / diff(get(gca, 'YLim')), ...
     'VerticalAlignment', 'bottom')
text(xt(end) / 5, yt(end) * 0.8, sprintf('\\alpha = %4.3f', p(1) / p(2)), 'FontSize', fs, 'EdgeColor', 'k', 'Margin', 2)
title('Alpha (\alpha) factor determination', 'FontSize', fs)
xlabel('1/R_I', 'FontSize', fs)
ylabel('1/(R_F+1)', 'FontSize', fs)
fname = sprintf('images/%03d_alphafactor.png', u); u = u + 1;
set(f1,'PaperPositionMode','auto')
print(f1, fname, '-dpng');

channels.alpha = p(1) / p(2);
channels.eaca_edcd = p(2);



%% Calculate donor intensity ID, accetor intensity IA and FRET efficiency E
% Donor intensity:

alpha = channels.alpha;
S4m = channels.S4med;
close all

%% Overlap of FRET efficiencies
% Create a figure
f1 = figure('Units', 'Pixels', 'Position', [100 100 600 600]);
% Create axes
axes('Units', 'Pixels', 'Position', [110 70 450 450], 'FontSize', fs);

leg = {'7aa', '19aa', '32aa'};
for i = 1 : numel(channels.transfer)
    F1 = channels.transfer(i).F1 - logOff;
    F1 = F1(channels.transfer(i).posGate);
    F1 = F1(channels.transfer(i).nonnegGate);
    F2 = channels.transfer(i).F2 - logOff;
    F2 = F2(channels.transfer(i).posGate);
    F2 = F2(channels.transfer(i).nonnegGate);
    F3 = channels.transfer(i).F3 - logOff;
    F3 = F3(channels.transfer(i).posGate);
    F3 = F3(channels.transfer(i).nonnegGate);

    channels.transfer(i).ID = ( -1 * alpha * ( S1m - S2m * S3m ) * (F1 * S2m - F2 * S4m) + S1m * S2m * (F1 * S1m + F3 * S2m - F1 * S2m * S3m - F3 * S1m * S4m + F2 * (-1 + S3m * S4m) ) ) / (alpha * (S1m - S2m * S3m) * (-1 * S2m + S1m * S4m) );
    channels.transfer(i).IA = ( F3 * S1m - F2 * S3m ) / ( S1m - S2m * S3m );
    channels.transfer(i).E = ( S1m * S2m * ( F1 * ( S1m - S2m * S3m ) + F3 * ( S2m - S1m * S4m ) + F2 * ( -1 + S3m * S4m ) ) ) ./ ( -1 * alpha * ( S1m - S2m * S3m ) * ( F1 * S2m - F2 * S4m ) + S1m * S2m * ( F1 * S1m + F3 * S2m - F1 * S2m * S3m - F3 * S1m * S4m + F2 * ( -1 + S3m * S4m ) ) );
    channels.transfer(i).E2 = ( S2m * (F2 - F1 * S1m - F3 * S2m ) ) ./ ( alpha * F1 * S2m + S2m * ( F2 - F1 * S1m - F3 * S2m ) );

    channels.transfer(i).ERI = channels.transfer(i).RI ./ (alpha + channels.transfer(i).RI);
    channels.transfer(i).ERF = (channels.transfer(i).RF - 1) * channels.eaca_edcd;
    % Create a figure
    figure('Units', 'Pixels', 'Position', [100 100 600 600])
    % Create axes
    axes('Units', 'Pixels', 'Position', [110 70 450 450], 'FontSize', fs); %#ok<LAXES>

    E = channels.transfer(i).E;
    in = E > -0.20 & E < 1 & F1 > 100 & F1 < 1000000;
    dscatter(100 * E(in), F1(in), 'logy', 1)
    set(gca, 'XLim', [-20, 100], 'YLim', [100, 5e4 ^ ceil(max(F1(in)) / 5e4)])
    title(sprintf('FRET Efficiency Scatter Plot for %s', leg{i}), 'FontSize', fs)
    xlabel('FRET efficiency [%]', 'FontSize', fs)
    ylabel('GFP Intensity', 'FontSize', fs)
    axis square
    box on
    % Convert Y-axis ticks to legible format
    ytick = get(gca, 'YTick');
    yticklabel = cell(size(ytick));
    for j = 1 : numel(ytick)
        [V, E] = num2eng(ytick(j), 3);
        yticklabel{j} = strcat(num2str(round(str2double(V))), E); 
    end
    set(gca, 'YTickLabel', yticklabel);
    fname = sprintf('images/%03d_FRETeff%sLog.png', u, leg{i}); u = u + 1;
    set(gcf,'PaperPositionMode','auto')
    print(gcf, fname, '-dpng');
    set(gca, 'YScale', 'linear');
    set(gca, 'YTickMode', 'auto')
    % Convert Y-axis ticks to legible format
    ytick = get(gca, 'YTick');
    yticklabel = cell(size(ytick));
    for j = 1 : numel(ytick)
        [V, E] = num2eng(ytick(j), 3);
        yticklabel{j} = strcat(num2str(round(str2double(V))), E); 
    end
    set(gca, 'YTickLabel', yticklabel);
    fname = sprintf('images/%03d_FRETeff%sLin.png', u, leg{i}); u = u + 1;
    print(gcf, fname, '-dpng');

    figure(f1);
    E = channels.transfer(i).ERI;
    dscatter(100 * E(in), F1(in), 'logy', 1)
    hold on
end

% Get the color of the scatter plots right
h = findobj(gca, 'type', 'hggroup');

map = vertcat(cool(64), autumn(64), winter(64));
m = size(map, 1) / numel(h); % 64-elements is each colormap

Z = cell(1, numel(h));
C = Z;
cmin = zeros(size(Z));
cmax = cmin;
for i = 1 : numel(Z)
    Z{i} = get(h(i), 'CData');
    cmin(i) = min(Z{i}(:));
    cmax(i) = max(Z{i}(:));
end

% CData for surface
for i = 1 : numel(Z)
    C{i} = min(m, 1 + round((m - 1) * (Z{i} - cmin(i)) / ...
               (cmax(i) - cmin(i)))) + (i - 1) * m;
end

% CData for pcolor
for i = 1 : numel(Z)
    set(h(i), 'CData', C{i});
end
caxis([min(C{1}(:)) max(C{i}(:))])

legend(leg, 'Location', 'NorthWest')
set(gca, 'XLim', [-20, 100], 'YLim', [100, 5e4 ^ ceil(max(F1(in)) / 5e4)])
title('Fully Corrected FRET Efficiency Scatter Plot', 'FontSize', fs)
xlabel('FRET efficiency [%]', 'FontSize', fs)
ylabel('GFP Intensity', 'FontSize', fs)
axis square
box on
% Convert Y-axis ticks to legible format
ytick = get(gca, 'YTick');
yticklabel = cell(size(ytick));
for j = 1 : numel(ytick)
    [V, E] = num2eng(ytick(j), 3);
    yticklabel{j} = strcat(num2str(round(str2double(V))), E); 
end
set(gca, 'YTickLabel', yticklabel);
fname = sprintf('images/%03d_FRETeffLog.png', u); u = u + 1;
set(gcf,'PaperPositionMode','auto')
print(gcf, fname, '-dpng');
set(gca, 'YScale', 'linear');
set(gca, 'YTickMode', 'auto')
% Convert Y-axis ticks to legible format
ytick = get(gca, 'YTick');
yticklabel = cell(size(ytick));
for j = 1 : numel(ytick)
    [V, E] = num2eng(ytick(j), 3);
    yticklabel{j} = strcat(num2str(round(str2double(V))), E); 
end
set(gca, 'YTickLabel', yticklabel);
fname = sprintf('images/%03d_FRETeffLin.png', u); u = u + 1;
print(gcf, fname, '-dpng');

% Create a figure
figure('Units', 'Pixels', 'Position', [100 100 600 600])
% Create axes
axes('Units', 'Pixels', 'Position', [110 70 450 450], 'FontSize', fs);
h(1) = plot(-20:100, histc(channels.transfer(1).ERI * 100, -20:100), 'k-', 'LineWidth', 2);
hold on
h(2) = plot(-20:100, histc(channels.transfer(2).ERI * 100, -20:100), 'k-.', 'LineWidth', 2);
h(3) = plot(-20:100, histc(channels.transfer(3).ERI * 100, -20:100), 'k--', 'LineWidth', 2);
title('R_I-Based FRET Efficiency', 'FontSize', fs)
xlabel('FRET Efficiency [%]', 'FontSize', fs)
ylabel('Frequency', 'FontSize', fs)
axis square
set(gca, 'XLim', [-20, 100])
box on
legend(h, leg, 'Location', 'NorthWest');
set(gcf,'PaperPositionMode','auto')
fname = sprintf('images/%03d_FRETeff.png', u); u = u + 1;
print(gcf, fname, '-dpng');

% Create a figure
figure('Units', 'Pixels', 'Position', [100 100 600 600])
% Create axes
axes('Units', 'Pixels', 'Position', [110 70 450 450], 'FontSize', fs);
h(1) = plot(-20:100, histc(channels.transfer(1).ERF * 100, -20:100), 'k-', 'LineWidth', 2);
hold on
h(2) = plot(-20:100, histc(channels.transfer(2).ERF * 100, -20:100), 'k-.', 'LineWidth', 2);
h(3) = plot(-20:100, histc(channels.transfer(3).ERF * 100, -20:100), 'k--', 'LineWidth', 2);
for i = 1 : 3
    ERF = sort(channels.transfer(i).ERF);
    ERF = 100 * ERF(round(numel(ERF) * 0.05) : round(numel(ERF) * 0.95));
    fprintf('ERF %s: %4.3f, %4.3f, %4.3f, %4.3f\n', leg{i}, mean(ERF), median(ERF), std(ERF), std(ERF) / sqrt(numel(ERF)));
    %plot(-20:100, histc(ERF, -20:100), ll{i})
    ERI = sort(channels.transfer(i).ERI);
    ERI = 100 * ERI(round(numel(ERI) * 0.05) : round(numel(ERI) * 0.95));
    fprintf('ERI %s: %4.3f, %4.3f, %4.3f, %4.3f\n', leg{i}, mean(ERI), median(ERI), std(ERI), std(ERI) / sqrt(numel(ERI)));
end
title('R_F-Based FRET Efficiency', 'FontSize', fs)
xlabel('FRET Efficiency [%]', 'FontSize', fs)
ylabel('Frequency', 'FontSize', fs)
axis square
set(gca, 'XLim', [-20, 100])
box on
legend(h, leg, 'Location', 'NorthWest');
set(gcf,'PaperPositionMode','auto')
fname = sprintf('images/%03d_FRETeffSim.png', u); u = u + 1;
print(gcf, fname, '-dpng');


%% Create a 5-color scatter plot
% Create a figure
figure('Units', 'Pixels', 'Position', [100 100 600 600])
% Create axes
axes('Units', 'Pixels', 'Position', [110 70 450 450], 'FontSize', fs);
hold on

x = logspace(1, 6, 200);
plots = {'donor', 'acceptor', 'transfer', 'transfer', 'transfer'; ...
         1, 1, 1, 2, 3};
hh = zeros(numel(x), numel(x), size(plots, 2));
for i = 1 : size(plots, 2);
    gate = channels.(plots{1, i})(plots{2, i}).posGate;
    sc = hist3([channels.(plots{1, i})(plots{2, i}).F1(gate), ...
                channels.(plots{1, i})(plots{2, i}).F2(gate)], {x, x});
    sc(sc < max(sc(:)) / 100) = NaN;
    hh(:, :, i) = sc;
end


[xx, yy] = meshgrid(x, x);

hs = zeros(1, size(plots, 2));
for i = 1 : numel(hs)
    hs(i) = surf(yy', xx', squeeze(hh(:, :, i))');
end

%h = findobj(gca, 'type', 'hggroup');
% colormap
map = vertcat(copper, autumn, winter, cool, bone);
colormap(map);
% Initially, both CDatas are equal to Z.
m = size(map, 1) / numel(hs); % 64-elements is each colormap

Z = cell(1, numel(hs));
C = Z;
cmin = zeros(size(Z));
cmax = cmin;
for i = 1 : numel(Z)
    Z{i} = get(hs(i), 'CData');
    cmin(i) = min(Z{i}(:));
    cmax(i) = max(Z{i}(:));
end

% CData for surface
for i = 1 : numel(Z)
    C{i} = min(m, 1 + round((m - 1) * (Z{i} - cmin(i)) / ...
               (cmax(i) - cmin(i)))) + (i - 1) * m;
end

% CData for pcolor
for i = 1 : numel(Z)
    set(hs(i), 'CData', C{i});
end
caxis([min(C{1}(:)) max(C{i}(:))])
set(gca, 'XLim', [min(x) max(x)], 'YLim', [min(x), max(x)], 'YScale', 'log', 'XScale', 'log')
set(hs, 'FaceAlpha', 0.1, 'EdgeColor', 'none')





figure('Units', 'Pixels', 'Position', [100 100 600 600])
u = 0;

for i = 1 : numel(channels.transfer)
    F1 = channels.transfer(i).F1 - logOff;
    F1 = F1(channels.transfer(i).posGate);
    F1 = F1(channels.transfer(i).nonnegGate);
    E = channels.transfer(i).E;
    in = E > -0.20 & E < 1 & F1 > 100 & F1 < 1000000;
    axes('Units', 'Pixels', 'Position', [110 520 - i * 150 450 150], 'FontSize', fs);
    dscatter(100 * E(in), F1(in), 'logy', 1)
    set(gca, 'XLim', [-20 100], 'YLim', [300 30000], 'Box', 'on', ...
        'XTick', -20 : 20 : 100, 'XTickLabel', {'' '' '' '' '' '' '' }, ...
        'YTick', [1000 10000], 'YTickLabel', {'1k', '10k'})

    text(95, 20000, leg{i}, 'FontSize', fs, 'HorizontalAlignment', 'right', ...
         'VerticalAlignment', 'top')
    if i == 2
        ylabel('GFP Intensity [Counts]', 'FontSize', fs);
    end
end

set(gca, 'XTickLabel', -20 : 20 : 100)
xlabel('FRET Efficiency [%]', 'FontSize', fs);

colormap('autumn')

set(gcf,'PaperPositionMode','auto')
fname = sprintf('images/%03d_FRETeffPan.png', u); u = u + 1;
print(gcf, fname, '-dpng');
saveas(gcf, [fname(1 : end - 3), 'fig'])

end

%% SCATTER_PLOT
function gate = scatter_plot(par)
% gate = scatter_plot(par)
%
%   Draw scatter plots and gate the data
%
%   gate:   Mask of the selected cells
%   
%   par.xdata:  X coordinate data set
%   par.ydata:  Y coordinate data set
%   par.xlim:   X axis limits (optional)
%   par.ylim:   Y axis limits (optional)
%   par.type:   Type of data
%               0:  background
%               1:  acceptor
%               2:  donor
%               3:  donor-acceptor
%   par.gate:   Gate coordinates (optional)
%   par.quadrant:   Quadrant values (optional)
%   par.mergeq:     Merge quadrants (optional)
%   par.title:  First word(s) of the title (optional)
%   par.colormap:   Colormap of the gated area (optional)
%   par.name:   Text in the title before the word Scatter (optional)
%   par.xaxis:  X axis title (optional)
%   par.yaxis:  Y axis title (optional)
%   par.dscatter:   dscatter parameters (optional)
%   par.fname:  input dataset filename (optional)
%   par.fiden:  filename identifier (optional)
%   par.findex: filename index (optional)

global fs    % Font Size

if ~isfield(par, 'name')
    par.name = '';
else
    % add a space if there is any content
    if ~isempty(par.name)
        par.name = horzcat(par.name, ' ');
    end
end

%close old figures
close all
% start a new figure
figure('Units', 'Pixels', 'Position', [100 100 600 600]);

% Create dscatter parameters if none are passed
if ~isfield(par, 'dscatter')
    par.dscatter = {};
end

% Create new axes
axes('Units', 'Pixels', 'Position', [110 70 450 450]);
% Create a scatter plot
dscatter(par.xdata, par.ydata, par.dscatter{:});

% Set a box around the plot
box on

% Set graph limits
if isfield(par, 'xlim')
    set(gca, 'Xlim', par.xlim)
end
if isfield(par, 'ylim')
    set(gca, 'Ylim', par.ylim)
end

% create a square plot
axis square

% Label the X-axis
if isfield(par, 'xaxis')
    xlabel(par.xaxis, 'FontSize', fs)
end
if isfield(par, 'yaxis')
    ylabel(par.yaxis, 'FontSize', fs)
end

% Set the font size
set(gca, 'FontSize', fs);

% Give the graph a title
switch par.type
    case 0
        par.name = horzcat(par.name, 'Background');
    case 1
        par.name = horzcat(par.name, 'Acceptor');
    case 2
        par.name = horzcat(par.name, 'Donor');
    case 3
        par.name = horzcat(par.name, 'Transfer');
end
if isfield(par, 'title')
    par.title = sprintf('%s Scatter Plot of %s Sample', ...
                        par.title, par.name);
else
    par.title = sprintf('Scatter Plot of %s Sample', par.name);
end
ht = title(par.title, 'FontSize', fs, 'Units', 'pixels');
% Break line in the title if too long
if min(get(ht, 'Extent')) < -10
    % find space nearest to the center
    in = regexp(par.title, ' ');
    [~, in1] = min(abs(in - numel(par.title) / 2));
    par.title = {par.title(1 : in(in1) - 1), par.title(in(in1) + 1 : end)};
	set(ht, 'String', par.title)
end


% Make sure the logarithmic scale has got all the points
if strcmpi(get(gca, 'XScale'), 'log')
    xl = get(gca, 'XLim');
    xl(1) = max(xl(1), 1);
    set(gca, 'XTick', 10 .^ (ceil(log10(xl(1))) : floor(log10(xl(2)))))
end
if strcmpi(get(gca, 'YScale'), 'log')
    yl = get(gca, 'YLim');
    yl(1) = max(yl(1), 1);
    set(gca, 'YTick', 10 .^ (ceil(log10(yl(1))) : floor(log10(yl(2)))))
end

% Convert X-axis ticks to legible format
xtick = get(gca, 'XTick');
xticklabel = cell(size(xtick));
for i = 1 : numel(xtick)
    [V, E] = num2eng(xtick(i), 3);
    xticklabel{i} = strcat(num2str(round(str2double(V))), E); 
end
set(gca, 'XTickLabel', xticklabel);

% Convert Y-axis ticks to legible format
ytick = get(gca, 'YTick');
yticklabel = cell(size(ytick));
for i = 1 : numel(ytick)
    [V, E] = num2eng(ytick(i), 3);
    yticklabel{i} = strcat(num2str(round(str2double(V))), E); 
end
set(gca, 'YTickLabel', yticklabel);

% label the graph with the filename in the lower right corner
if isfield(par, 'fname')
    % create an empty invisible window
    ha = gca;
    xl = get(gca, 'XLim');
    yl = get(gca, 'YLim');
    ha1 = axes;
    axis square
    set(ha1, 'Units', 'Pixels', 'Position', get(ha, 'Position'), ...
        'Color', 'none', 'XLim', xl, 'YLim', yl, 'XTick', [], 'YTick', [])
    text(xl(1) + diff(xl) * 0.98, yl(1) + diff(yl) * 0.02, ...
         regexprep(par.fname, '_', '\\_'), ...
         'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
         'FontSize', fs * 0.7);
	axes(ha)
end

% Set a colormap for the graph
map = copper(128);
map = map(65 : end, :);
colormap(map)

% Create the gate
gate = true(size(par.xdata));

% Plot a gate and adjust graph colormap if a gate is set
if isfield(par, 'gate')
    % Create an enclosed gate polygon "G"
    G = vertcat(par.gate, par.gate(1, :));

    % plot the gate "G"
    hold on
    if isfield(par, 'colormap')
        col = eval(par.colormap);
        col = col(1, :);
    else
        col = 'b';
    end
    plot(G(:, 1), G(:, 2), '-', 'LineWidth', 2, 'Color', col);

    % select cells within the gate "G"
    for i = 1 : (size(G, 1) - 1)
        gate = gate & ...
               (G(i + 1, 1) - G(i, 1)) * (par.ydata - G(i, 2)) - ...
               (G(i + 1, 2) - G(i, 2)) * (par.xdata - G(i, 1)) > 0;
    end

    if any(gate)
        % Create a scatter plot of the gated cells
        dscatter(par.xdata(gate), par.ydata(gate), par.dscatter{:});

        % Get the color of the scatter plots right
        h = findobj(gca, 'type', 'hggroup');

        if isfield(par, 'colormap')
            map = vertcat(feval(par.colormap, 64), map);
        else
            map = vertcat(winter(64), map);
        end
        colormap(map)

        % Initially, both CDatas are equal to Z.
        m = size(map, 1) / 2; % 64-elements is each colormap

        Z1 = get(h(1), 'CData'); c1min = min(Z1(:)); c1max = max(Z1(:));
        Z2 = get(h(2), 'CData'); c2min = min(Z2(:)); c2max = max(Z2(:));

        % CData for surface
        C1 = min(m,round((m-1)*(Z1-c1min)/(c1max-c1min))+1);
        C2 = min(m,round((m-1)*(Z2-c2min)/(c2max-c2min))+1) + m;
        % CData for pcolor
        set(h(1),'CData',C1);
        set(h(2),'CData',C2);
        caxis([min(C1(:)) max(C2(:))])
    end
end

% Plot quadrants and adjust graph colormap if a gate is set
if isfield(par, 'quadrant')
    % Create the quadrant
    gate = true(size(par.xdata, 1), 4);

    % plot the quadrants
    hold on

    % Graph axis limits
    xl = get(gca, 'XLim');
    if strcmpi(get(gca, 'XScale'), 'log')
        xl(1) = max(xl(1), 1);
    end
    yl = get(gca, 'YLim');
    if strcmpi(get(gca, 'YScale'), 'log')
        yl(1) = max(yl(1), 1);
    end
    plot(par.quadrant(1) * [1 1], [yl(1), par.quadrant(2)], 'k-', 'LineWidth', 2);
    plot(par.quadrant(1) * [1 1], [par.quadrant(2), yl(2)], 'k-', 'LineWidth', 2);
    plot([xl(1), par.quadrant(1)], par.quadrant(2) * [1 1], 'k-', 'LineWidth', 2);
    plot([par.quadrant(1), xl(2)], par.quadrant(2) * [1 1], 'k-', 'LineWidth', 2);

    % select cells within the quadrants
    gate(:, 1) = ...
        par.xdata <= par.quadrant(1) & par.ydata <= par.quadrant(2);
    gate(:, 2) = ...
        par.xdata > par.quadrant(1) & par.ydata <= par.quadrant(2);
    gate(:, 3) = ...
        par.xdata <= par.quadrant(1) & par.ydata > par.quadrant(2);
    gate(:, 4) = ...
        par.xdata > par.quadrant(1) & par.ydata > par.quadrant(2);

    % Check if there are any quadrants to merge
    if ~isfield(par, 'mergeq')
        par.mergeq = 1 : 4;     % each quadrant is unique, on its own
    end

    % At least two cells in a gate
    gated = [true false false false];
    for i = unique(par.mergeq(2 : 4))
        G = sum(gate(:, par.mergeq == i), 2) > 0;
        if sum(G) > 1
            % Create a scatter plot of the gated cells
            dscatter(par.xdata(G), par.ydata(G), par.dscatter{:});
            gated(find(par.mergeq == i, 1, 'first')) = true;
        end
    end

    % Get the color of the scatter plots right
    h = findobj(gca, 'type', 'hggroup');

    if isfield(par, 'colormap')
        if iscell(par.colormap)
            for k = 1 : 3
                if ischar(par.colormap{k})
                    par.colormap{k} = feval(par.colormap{k}, 64);
                end
            end
            map = vertcat(par.colormap{gated(2 : 4)}, map);
        else
            if ischar(par.colormap)
                par.colormap = ...
                    repmat(feval(par.colormap, 64), sum(gated(2 : 4)), 1);
            end
            map = vertcat(par.colormap, map);
        end
    else
        par.colormap = {'winter', 'summer', 'autumn'};
        par.colormap = cellfun(@(x) feval(x, 64), par.colormap, 'UniformOutput', false);
        map = vertcat(par.colormap{gated(2 : 4)}, map);
    end
    colormap(map)

    % Initially, both CDatas are equal to Z.

    m = size(map, 1) / numel(h); % 64-elements is each colormap

    Z = cell(1, numel(h));
    C = Z;
    cmin = zeros(size(Z));
    cmax = cmin;
    for i = 1 : numel(Z)
        Z{i} = get(h(i), 'CData');
        cmin(i) = min(Z{i}(:));
        cmax(i) = max(Z{i}(:));
    end

    % CData for surface
    for i = 1 : numel(Z)
        C{i} = min(m, 1 + round((m - 1) * (Z{i} - cmin(i)) / ...
                   (cmax(i) - cmin(i)))) + (i - 1) * m;
    end

    % CData for pcolor
    for i = 1 : numel(Z)
        set(h(i), 'CData', C{i});
    end
    caxis([min(C{1}(:)) max(C{i}(:))])
    set(gca, 'XLim', xl, 'YLim', yl)
    
    % Get positions of corners to write into
    pos = vertcat(diff(xl) * [0.02, 0.98, 0.02, 0.98] + xl(1), ...
                  diff(yl) * [0.02, 0.02, 0.98, 0.98] + yl(1));
    if isfield(par, 'fname')
        axes(ha1)
        pos(2, 2) = diff(yl) * 0.06 + yl(1);
    else
        ha = gca;
        axes;
        axis square
        set(ha1, 'Position', get(ha, 'Position'), 'Color', 'none', ...
            'XLim', xl, 'YLim', yl, 'XTick', [], 'YTick', [])
    end

    % Printe percentages of cells in gates
    align = {'left', 'right', 'left', 'right'; ...
             'bottom', 'bottom', 'top', 'top'};
    for i = 1 : 4
        perc = sum(gate(:, i)) / size(gate, 1) * 100;
        if perc == 100;
            text(pos(1, i), pos(2, i), ...
                sprintf('100 %% (%d)', sum(gate(:, i))), ...
                'HorizontalAlignment', align{1, i}, ...
                'VerticalAlignment', align{2, i}, 'FontSize', 0.7 * fs);
        else
            text(pos(1, i), pos(2, i), ...
                sprintf('%.02f %% (%d)', perc, sum(gate(:, i))), ...
                'HorizontalAlignment', align{1, i}, ...
                'VerticalAlignment', align{2, i}, 'FontSize', 0.7 * fs);
        end
    end
end

if exist('ha1', 'var')
    if ishandle(ha1)
        axes(ha1)
    end
end


% Save the graph
if isfield(par, 'fiden')
    par.name = regexprep(par.name, ' ', '-');
    par.fiden = regexprep(par.fiden, ' ', '');
    if isfield(par, 'findex')
        fname = sprintf('images%s%03d_%s_%s.png', ...
                        filesep, par.findex, par.fiden, par.name);
    else
        fname = sprintf('images%s%s_%s.png', filesep, par.fiden, par.name);
    end
    set(gcf,'PaperPositionMode','auto')
    print(gcf, fname, '-dpng');
end

% end the function "gate = scatter_plot(par)"
end

%% HISTOGRAM
function gate = histogram(par)
% gate = histogram(par)
%
%   Draw histogram and gate the data
%
%   gate:   Mask of the selected cells
%   
%   par.data:   Input data set
%   par.xbins:  Bin vector or numbe of bins
%   par.xlim:   X axis limits (optional)
%   par.ylim:   Y axis limits (optional)
%   par.type:   Type of data
%               0:  background
%               1:  acceptor
%               2:  donor
%               3:  donor-acceptor
%   par.gate:   Gate coordinates (optional)
%   par.title:  First word(s) of the title (optional)
%   par.color:  Color of the gated area (optional)
%   par.name:   Text in the title before the word Scatter (optional)
%   par.xaxis:  X axis title (optional)
%   par.yaxis:  Y axis title (optional)
%   par.fname:  input dataset filename (optional)
%   par.fiden:  filename identifier (optional)
%   par.findex: filename index (optional)
global fs

if ~isfield(par, 'name')
    par.name = '';
else
    % add a space if there is any content
    if ~isempty(par.name)
        par.name = horzcat(par.name, ' ');
    end
end

%close old figures
close all
% start a new figure
figure('Units', 'Pixels', 'Position', [100 100 600 600]);

% Create bins of the histogram
[H, par.xbins] = hist(par.data, par.xbins);

% Create a gate for the data
gate = true(size(par.data));

% Create single- or dual-color histogram
if isfield(par, 'gate')
    for i = 1 : size(par.gate, 1)
        gate = ...
            gate & par.data > par.gate(i, 1) & par.data < par.gate(i, 2);
    end
    % Plot a two-color bargraph
    H = vertcat(hist(par.data(~gate), par.xbins), ...
                hist(par.data(gate), par.xbins));
else
    H = vertcat(H, zeros(size(H)));
end

% Create axes
axes('Units', 'Pixels', 'Position', [110 70 450 450], 'FontSize', fs);

% Plot the histogram
hb = bar(par.xbins, H', 'stacked');
set(hb, 'BarWidth', 1, 'EdgeColor', 'none');

% Set the colors of the histogram
if isfield(par, 'color')
    col = par.color;
else
    col = 'r';
end
set(hb(1), 'FaceColor', [0 0 0.5]);
set(hb(2), 'FaceColor', col);

% Set a box around the plot
box on

% Set graph limits
if isfield(par, 'xlim')
    set(gca, 'Xlim', par.xlim)
end
if isfield(par, 'ylim')
    set(gca, 'Ylim', par.ylim)
end

% create a square plot
axis square

% Label the X-axis
if isfield(par, 'xaxis')
    xlabel(par.xaxis, 'FontSize', fs)
end
if isfield(par, 'yaxis')
    ylabel(par.yaxis, 'FontSize', fs)
end

% Give the graph a title
switch par.type
    case 0
        par.name = horzcat(par.name, 'Background');
    case 1
        par.name = horzcat(par.name, 'Acceptor');
    case 2
        par.name = horzcat(par.name, 'Donor');
    case 3
        par.name = horzcat(par.name, 'Transfer');
end
if isfield(par, 'title')
    par.title = sprintf('%s Scatter Plot of %s Sample', ...
                        par.title, par.name);
else
    par.title = sprintf('Scatter Plot of %s Sample', par.name);
end
% Give the plot a title
ht = title(par.title, 'FontSize', fs, 'Units', 'pixels');
% Break line in the title if too long
if min(get(ht, 'Extent')) < -10
    % find space nearest to the center
    in = regexp(par.title, ' ');
    [~, in1] = min(abs(in - numel(par.title) / 2));
    par.title = {par.title(1 : in(in1) - 1), par.title(in(in1) + 1 : end)};
	set(ht, 'String', par.title)
end


% % Convert X-axis ticks to legible format
% xtick = get(gca, 'XTick');
% xticklabel = cell(size(xtick));
% for i = 1 : numel(xtick)
%     [V, E] = num2eng(xtick(i), 3);
%     xticklabel{i} = strcat(num2str(round(str2double(V))), E); 
% end
% set(gca, 'XTickLabel', xticklabel);

% Convert Y-axis ticks to legible format
% ytick = get(gca, 'YTick');
% yticklabel = cell(size(ytick));
% for i = 1 : numel(ytick)
%     [V, E] = num2eng(ytick(i), 3);
%     yticklabel{i} = strcat(num2str(round(str2double(V))), E); 
% end
% set(gca, 'YTickLabel', yticklabel);

% label the graph with the filename in the lower right corner
if isfield(par, 'fname')
    xl = get(gca, 'XLim');
    yl = get(gca, 'YLim');
    text(xl(1) + diff(xl) * 0.98, yl(1) + diff(yl) * 0.98, ...
         regexprep(par.fname, '_', '\\_'), 'FontSize', 0.7 * fs, ...
         'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
end

% plot gates
if isfield(par, 'gate')
    hold on
    % plot the gate
    yl = get(gca, 'YLim');      % Y limits of the graph
    for i = 1 : size(par.gate, 1)
        x = par.gate(i, [1 1 2 2]);
        y = yl(2) * [0.97 0.99 0.99 0.97];
        plot(x, y, '-', 'LineWidth', 2, 'Color', col);
    end
end

% Save the graph
if isfield(par, 'fiden')
    par.name = regexprep(par.name, ' ', '-');
    if isfield(par, 'findex')
        fname = sprintf('images%s%03d_%s_%s.png', ...
                        filesep, par.findex, par.fiden, par.name);
    else
        fname = sprintf('images%s%s_%s.png', filesep, par.fiden, par.name);
    end
    set(gcf,'PaperPositionMode','auto')
    print(gcf, fname, '-dpng');
end

% end the function "gate = scatter_plot(par)"
end

