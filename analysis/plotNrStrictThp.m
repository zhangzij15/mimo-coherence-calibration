function files = plotNrStrictThp(csvPath,outputDir)
%PLOTNRSTRICTTHP Plot the short-SRS-age NR throughput confirmation.
%   FILES = PLOTNRSTRICTTHP(CSVPATH,OUTPUTDIR) reads confirm_thp.csv and
%   exports a vector PDF and a 300 dpi PNG from the same four-series figure.
arguments
    csvPath (1,1) string
    outputDir (1,1) string
end

required = {'snr_dB','ideal_thp_mbps','uncal_thp_mbps','feedback_thp_mbps', ...
    'reselect_thp_mbps','mc','selected_slots_per_mc','duty'};
data = readtable(csvPath);
assert(height(data)>=2 && all(ismember(required,data.Properties.VariableNames)), ...
    'UL256:THPPlotSchema','The confirmation table must contain the four throughput series and run metadata.');
for k = 1:numel(required)
    column = data.(required{k});
    assert(isnumeric(column) && isreal(column) && iscolumn(column) && all(isfinite(column)), ...
        'UL256:THPPlotSchema','Column %s must contain finite real numbers.',required{k});
end
assert(all(diff(data.snr_dB)>0),'UL256:THPPlotSNR', ...
    'Confirmation SNR values must be strictly increasing and unique.');
assert(all(data.mc>=20 & data.mc==fix(data.mc)) && all(data.mc==data.mc(1)), ...
    'UL256:THPPlotIncomplete','At least 20 complete paired MC realizations are required at each SNR.');
assert(all(data.selected_slots_per_mc>0 & data.selected_slots_per_mc==fix(data.selected_slots_per_mc)) && ...
    all(data.selected_slots_per_mc==data.selected_slots_per_mc(1)) && ...
    all(data.duty>0 & data.duty<=1) && all(abs(data.duty-data.duty(1))<1e-12), ...
    'UL256:THPPlotSchedule','The short-age slots and original-frame duty must agree across SNR points.');
thp = data{:,required(2:5)};
assert(all(thp>=0,'all'),'UL256:THPPlotSchema','Throughput must be nonnegative.');

if ~isfolder(outputDir), mkdir(outputDir); end
fig = figure('Visible','off','Color','w','Units','pixels','Position',[100 100 680 470]);
cleanup = onCleanup(@() close(fig));
layout = tiledlayout(fig,1,1,'TileSpacing','compact','Padding','loose');
ax = nexttile(layout);
hold(ax,'on');
colors = [0.13 0.15 0.17; 0.42 0.42 0.42; 0.00 0.38 0.66; 0.75 0.28 0.00];
lineStyles = {'-','--','-.',':'};
markers = {'^','v','o','s'};
curves = gobjects(1,4);
for k = 1:4
    curves(k) = plot(ax,data.snr_dB,thp(:,k), ...
        'Color',colors(k,:),'LineStyle',lineStyles{k}, ...
        'Marker',markers{k},'MarkerFaceColor','w','MarkerSize',5.5, ...
        'LineWidth',1.5);
end
set(ax,'FontName','Arial','FontSize',10,'Box','off','TickDir','out', ...
    'LineWidth',.8,'YGrid','on','GridAlpha',.16,'Layer','top');
xticks(ax,data.snr_dB);
margin = .04*diff(data.snr_dB([1 end]));
xlim(ax,[data.snr_dB(1)-margin data.snr_dB(end)+margin]);
xlabel(ax,'SNR (dB)');
ylabel(ax,'THP (Mbps)');
title(layout,'4T256R NR throughput for short SRS age', ...
    'FontName','Arial','FontSize',13,'FontWeight','bold');
lg = legend(ax,curves,{'Ideal coherent','Non-ideal', ...
    'Non-ideal + Feedback calibration','Non-ideal + Scheduling calibration'}, ...
    'NumColumns',2,'Orientation','horizontal','Box','off','FontSize',9);
lg.Layout.Tile = 'south';
xlabel(layout,{sprintf('SRS age ≤20 ms. MC %d with %d selected slots per realization.', ...
    data.mc(1),data.selected_slots_per_mc(1)), ...
    'THP uses the full 50 ms simulation interval. Scheduling calibration uses Reselect.'}, ...
    'FontName','Arial','FontSize',8.5,'Interpreter','none');
drawnow;
files = fullfile(outputDir,["nr_strict_thp_vs_snr.pdf";"nr_strict_thp_vs_snr.png"]);
exportgraphics(fig,files(1),'ContentType','vector','BackgroundColor','white');
exportgraphics(fig,files(2),'Resolution',300,'BackgroundColor','white');
end
