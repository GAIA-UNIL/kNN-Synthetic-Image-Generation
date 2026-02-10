function visualiseMetrics(nbImages,pixelWise,targetVar,climateVars,targetDim,refValidation,synImages,validationMetric,sortedDates,climateData,metricV,nanValue,varLegend,varRange,errRange,metricKNN,LdateStart,LdateEnd,QdateStart,QdateEnd,daysRange,outputTime,bootstrap,outDir,createGIF)

%
%
%
% REDO DOCUMENTATION
%
%
%

startLdate = datetime(LdateStart,'ConvertFrom','yyyyMMdd','Format','dd-MM-yyyy');
startQdate = datetime(QdateStart,'ConvertFrom','yyyyMMdd','Format','dd-MM-yyyy');
endLdate   = datetime(LdateEnd,'ConvertFrom','yyyyMMdd','Format','dd-MM-yyyy');
endQdate   = datetime(QdateEnd,'ConvertFrom','yyyyMMdd','Format','dd-MM-yyyy');

numPeriods = length(startQdate);

targetVarL = lower(targetVar);

for k = 1:numel(targetVar)
    if bootstrap == true
        dates = datetime(cell2mat(validationMetric.(targetVarL(k))(:,1)),'ConvertFrom','yyyyMMdd','Format','dd-MM-yyyy');
        % Validation data
        singleDataVal = cell2mat(validationMetric.(targetVarL(k))(:,3));
        maxValuesVal  = zeros(size(validationMetric.(targetVarL(k)), 1), 1);
        meanValuesVal = zeros(size(validationMetric.(targetVarL(k)), 1), 1);
        minValuesVal  = zeros(size(validationMetric.(targetVarL(k)), 1), 1);
        for i = 1:size(validationMetric.(targetVarL(k)), 1)
            valuesVal = sort(validationMetric.(targetVarL(k)){i, 2});  % Extract values from the second column of the cell
            maxValuesVal(i)  = valuesVal(end);  % Compute the maximum value
            meanValuesVal(i) = mean(valuesVal);
            minValuesVal(i)  = valuesVal(1);  % Compute the minimum value
        end
        % Synthetic data
        varBs = strcat(targetVar(k),'_stochastic');
        varBD = strcat(targetVar(k),'_Distances');
        refData    = squeeze(mean(mean(refValidation.(targetVarL(k)),1,'omitnan'),2,'omitnan'));  % Extract mean of ref variable
        synData    = squeeze(mean(mean(synImages.(targetVarL(k)),1,'omitnan'),2,'omitnan'));
        maxValues  = zeros(size(validationMetric.(targetVarL(k)), 1), 1);
        meanValues = zeros(size(validationMetric.(targetVarL(k)), 1), 1);
        minValues  = zeros(size(validationMetric.(targetVarL(k)), 1), 1);
        for i = 1:size(validationMetric.(targetVarL(k)), 1)
            synValues     = sort(squeeze(mean(mean(synImages.(varBs(k)){i},'omitnan'),'omitnan')));  % Extract mean of variable and sort
            maxValues(i)  = synValues(end);  % Compute the maximum value
            meanValues(i) = mean(synValues);
            minValues(i)  = synValues(1);  % Compute the minimum value
        end
        bestDist   = cellfun(@(x) min(x(:)), synImages.(varBD));

        % -----------------------------------------------------------------------------------------------------

        % === PLOT 1: Mean ===
        if numPeriods > 1
            for period = 1:numPeriods
                % Filter dates for current validation period
                periodIndices = dates >= startQdate(period) & dates <= endQdate(period);

                % Subset data for that period
                datesP = dates(periodIndices);
                synDataP = synData(periodIndices);
                refDataP = refData(periodIndices);
                maxValuesP = maxValues(periodIndices);
                minValuesP = minValues(periodIndices);

                % --- Plot ---
                figure();
                hold on
                inBetweenRegionX = [datesP', fliplr(datesP')];
                inBetweenRegionY = [maxValuesP', fliplr(minValuesP')];
                patch(inBetweenRegionX, inBetweenRegionY, 'k', 'LineStyle', 'none', 'FaceAlpha', 0.25)
                plot(datesP, synDataP, 'k-', 'LineWidth', 1)
                plot(datesP, refDataP, 'r-', 'LineWidth', 1)
                hold off

                % Compute metrics
                r = corr(synDataP, refDataP);
                alpha = std(synDataP) / std(refDataP);
                beta  = mean(synDataP) / mean(refDataP);
                kgeSynRef = 1 - sqrt((r-1)^2 + (alpha-1)^2 + (beta-1)^2);
                str = {['KGE: ' num2str(kgeSynRef,'%.5f')], ...
                    ['r: ' num2str(r,'%.5f') ', \alpha: ' num2str(alpha,'%.5f') ', \beta: ' num2str(beta,'%.5f')]};

                % Titles and labels
                title([convertStringsToChars(targetVar(k)) ' - MEAN'])
                subtitle(str)
                xlabel('Date')
                ylabel(['Mean ' convertStringsToChars(targetVar(k))])
                legend('Synthetic data spread','Deterministic mean','Reference data mean')
                grid on
                set(gcf, 'Color', 'white');

                % Save figure
                saveas(gcf, fullfile(outDir, ...
                    ['bsValidation_AVG_' convertStringsToChars(targetVar(k)) '_' char(startQdate(period)) '_' char(endQdate(period)) '.png']));
                close(gcf);
            end
        end

        % === PLOT 2: RMSE / metric-based ===
        if numPeriods > 1
            for period = 1:numPeriods
                % Filter dates for current validation period
                periodIndices = dates >= startQdate(period) & dates <= endQdate(period);

                % Subset data for that period
                datesP = dates(periodIndices);
                singleDataValP = singleDataVal(periodIndices);
                bestDistP = bestDist(periodIndices);
                maxValuesValP = maxValuesVal(periodIndices);
                minValuesValP = minValuesVal(periodIndices);

                % --- Plot ---
                figure();
                hold on
                inBetweenRegionX = [datesP', fliplr(datesP')];
                inBetweenRegionY = [maxValuesValP', fliplr(minValuesValP')];
                patch(inBetweenRegionX, inBetweenRegionY, 'k', 'LineStyle', 'none', 'FaceAlpha', 0.25)
                plot(datesP, singleDataValP, 'k-', 'LineWidth', 1)
                yyaxis right
                plot(datesP, bestDistP, 'r', 'LineWidth', 1)

                % Right Y label
                switch metricKNN
                    case 1, ylabel('RMSE')
                    case 2, ylabel('MAE')
                    case 3, ylabel('1-bSPEM')
                    case 4, ylabel('Hellinger distance')
                    case 5, ylabel('0.5*(1-bSPEM) + 0.5*Hellinger')
                    case 6, ylabel('SPAEF')
                end

                yyaxis left
                switch metricV
                    case 1, title([convertStringsToChars(targetVar(k)) ' - MAE'])
                    case 2, title([convertStringsToChars(targetVar(k)) ' - RMSE'])
                    case 3, title([convertStringsToChars(targetVar(k)) ' - SPEM'])
                    case 4, title([convertStringsToChars(targetVar(k)) ' - SPAEF'])
                    case 5, title([convertStringsToChars(targetVar(k)) ' - KGE'])
                end

                % Compute statistics
                str = {['Mean RMSE: ' num2str(mean(singleDataValP),'%.5f')], ...
                    ['RMSE - MAE correlation: ' num2str(corr(singleDataValP,bestDistP),'%.5f')]};
                subtitle(str)

                xlabel('Date')
                switch metricV
                    case 1, ylabel('MAE')
                    case 2, ylabel('RMSE')
                    case 3, ylabel('SPEM')
                    case 4, ylabel('SPAEF')
                    case 5, ylabel('KGE')
                end
                ax = gca;
                ax.YAxis(1).Color = 'k';
                ax.YAxis(2).Color = 'r';
                linkprop(ax.YAxis, 'Limits');
                legend('bootstrap ensembles','Deterministic','Best distance')
                grid on
                set(gcf, 'Color', 'white');
                
                yyaxis left
                ylim([0, max(ylim)])
                yyaxis right
                ylim([0, max(ylim)])

                % Save figure
                saveas(gcf, fullfile(outDir, ...
                    ['bsValidation_RMSE_' convertStringsToChars(targetVar(k)) '_' char(startQdate(period)) '_' char(endQdate(period)) '.png']));
                close(gcf);
            end
        end
    else
        dates = datetime(validationMetric.(targetVarL(k))(:,1),'ConvertFrom','yyyyMMdd','Format','dd-MM-yyyy');
        if targetDim == 1
            figure();
            plot(dates, refValidation.(targetVarL(k)),'Color','r','DisplayName','Reference');
            hold on
            plot(dates, synImages.(targetVarL(k)),'Color','b','DisplayName','Synthetic')
            hold off
            if length(startQdate) == 1
                if strcmp(startLdate, startQdate)
                    subtitle([['Learning period: ' endQdate '-' endLdate] str])
                elseif strcmp(endQdate, endLdate)
                    subtitle([['Learning period: ' startLdate '-' startQdate] str])
                else
                    subtitle([['Learning period: ' startLdate '-' startQdate ' - ' endQdate '-' endLdate] str])
                end
            end
            xlabel('Date')
            ylabel('Discharge [m^{3}/s]')
            title(targetVar(k))
            legend()
            set(gcf, 'color', 'white');
            grid on
            saveas(gcf,strcat(outDir,['\' convertStringsToChars(targetVar(k)) '.png']))
        end

        % --------------------------------------------------------------------
        % === Plot for each individual validation period ===
        if numPeriods > 1
            for period = 1:numPeriods
                figure();
                % Filter data for current validation period
                periodIndices = dates >= startQdate(period) & dates <= endQdate(period);
                periodData = validationMetric.(targetVarL(k))(periodIndices, :);
                % Plotting
                plot(datetime(periodData(:,1), 'ConvertFrom', 'yyyyMMdd', 'Format', 'dd/MM/yyyy'), periodData(:,2));
                yline(mean(periodData(:,2)), '-', ['Mean: ' num2str(mean(periodData(:,2)))], 'Color', 'r');
                % Title
                switch metricV
                    case 1, title([convertStringsToChars(targetVar(k)) ' - MAE (mean: ' num2str(mean(periodData(:,2))) ')']);
                    case 2, title([convertStringsToChars(targetVar(k)) ' - RMSE (mean: ' num2str(mean(periodData(:,2))) ')']);
                    case 3, title([convertStringsToChars(targetVar(k)) ' - SPEM (mean: ' num2str(mean(periodData(:,2))) ')']);
                    case 4, title([convertStringsToChars(targetVar(k)) ' - SPAEF (mean: ' num2str(mean(periodData(:,2))) ')']);
                    case 5, title([convertStringsToChars(targetVar(k)) ' - KGE (mean: ' num2str(mean(periodData(:,2))) ')']);
                end
                xlabel('Date');
                switch metricV
                    case 1, ylabel('MAE');
                    case 2, ylabel('RMSE');
                    case 3, ylabel('SPEM');
                    case 4, ylabel('SPAEF');
                    case 5, ylabel('KGE');
                end
                set(gcf, 'color', 'white');
                grid on;
                % Save figure
                saveas(gcf, fullfile(outDir, ...
                    ['validation_' convertStringsToChars(targetVar(k)) '_' char(startQdate(period)) '_' char(endQdate(period)) '.png']));
                close(gcf);
            end
        end
        % === Plot for the total validation period (from first start to last end) ===
        figure();
        % Compute overall period range
        totalStart = startQdate(1);
        totalEnd   = endQdate(end);
        % Filter data for total period
        totalIndices = dates >= totalStart & dates <= totalEnd;
        totalData = validationMetric.(targetVarL(k))(totalIndices, :);
        % Plotting
        plot(datetime(totalData(:,1), 'ConvertFrom', 'yyyyMMdd', 'Format', 'dd/MM/yyyy'), totalData(:,2));
        yline(mean(totalData(:,2)), '-', ['Mean: ' num2str(mean(totalData(:,2)))], 'Color', 'r');
        % Title
        switch metricV
            case 1, title([convertStringsToChars(targetVar(k)) ' - MAE (mean: ' num2str(mean(totalData(:,2))) ')']);
            case 2, title([convertStringsToChars(targetVar(k)) ' - RMSE (mean: ' num2str(mean(totalData(:,2))) ')']);
            case 3, title([convertStringsToChars(targetVar(k)) ' - SPEM (mean: ' num2str(mean(totalData(:,2))) ')']);
            case 4, title([convertStringsToChars(targetVar(k)) ' - SPAEF (mean: ' num2str(mean(totalData(:,2))) ')']);
            case 5, title([convertStringsToChars(targetVar(k)) ' - KGE (mean: ' num2str(mean(totalData(:,2))) ')']);
        end
        xlabel('Date');
        switch metricV
            case 1, ylabel('MAE');
            case 2, ylabel('RMSE');
            case 3, ylabel('SPEM');
            case 4, ylabel('SPAEF');
            case 5, ylabel('KGE');
        end
        set(gcf, 'color', 'white');
        grid on;
        % Save total period plot
        saveas(gcf, fullfile(outDir, ...
            ['validation_' convertStringsToChars(targetVar(k)) '_total_' char(totalStart) '_' char(totalEnd) '.png']));
        close(gcf);

        % --------------------------------------------------------------------
        if targetDim ~= 1
            refData = refValidation.(targetVarL(k));
            synData = synImages.(targetVarL(k));
            % Handle nan values
            if ~isnan(nanValue)
                refData(refData == nanValue) = nan;
                synData(isnan(refData)) = nan;  % Mask synthetic data where reference is NaN
            end
            % === Plot for each individual validation period ===
            if numPeriods > 1
                for period = 1:numPeriods
                    % Time filtering
                    periodIndices = dates >= startQdate(period) & dates <= endQdate(period);
                    % Compute spatial means over time for current period
                    meanRefData = squeeze(mean(mean(refData(:,:,periodIndices), 1, 'omitnan'), 2, 'omitnan'));
                    meanSynData = squeeze(mean(mean(synData(:,:,periodIndices), 1, 'omitnan'), 2, 'omitnan'));
                    % Plot
                    figure();
                    plot(dates(periodIndices), meanRefData, 'r-', ...
                        dates(periodIndices), meanSynData, 'k-');
                    legend('Reference', 'Synthetic', 'Location', 'southeast');
                    xlabel('Date');
                    ylabel(varLegend);
                    title(['Mean ' convertStringsToChars(targetVarL(k))]);
                    % KGE and stats
                    r = corr(meanSynData, meanRefData, 'rows', 'complete');
                    alpha = std(meanSynData, 'omitnan') / std(meanRefData, 'omitnan');
                    beta = mean(meanSynData, 'omitnan') / mean(meanRefData, 'omitnan');
                    kgeSynRef = 1 - sqrt((r - 1)^2 + (alpha - 1)^2 + (beta - 1)^2);
                    str = {['KGE: ' num2str(kgeSynRef, '%.5f')], ...
                        ['r: ' num2str(r, '%.5f') ', \alpha: ' num2str(alpha, '%.5f') ', \beta: ' num2str(beta, '%.5f')]};
                    subtitle(str);
                    grid on; box off;
                    set(gcf, 'color', 'white');
                    % Save figure
                    saveas(gcf, fullfile(outDir, ...
                        ['correlation_' convertStringsToChars(targetVar(k)) '_' char(startQdate(period)) '_' char(endQdate(period)) '.png']));
                    close(gcf);
                end
            end
            % === Plot for the total period (overlay by year using DOY) ===
            totalStart = startQdate(1);
            totalEnd = endQdate(end);
            totalIndices = dates >= totalStart & dates <= totalEnd;
            % Compute spatial means over time for full period
            meanRefData = squeeze(mean(mean(refData(:,:,totalIndices), 1, 'omitnan'), 2, 'omitnan'));
            meanSynData = squeeze(mean(mean(synData(:,:,totalIndices), 1, 'omitnan'), 2, 'omitnan'));
            selDates = dates(totalIndices);
            % Extract year and day-of-year (DOY)
            [yearVec, ~, ~] = ymd(selDates);
            doy = day(selDates, 'dayofyear');
            % Unique years
            uniqueYears = unique(yearVec);
            numYears = numel(uniqueYears);
            figure('Units','inches','Position',[1 1 6 3.5])
            hold on;
            % Preallocate handles for legend (one handle per year — use synthetic line handle)
            plottedIdx = 1:2:numYears;   % store which years were actually plotted
            h = gobjects(length(plottedIdx),1);
            % Prepare color map
            cmap = turbo(length(plottedIdx)); % distinct colors per year
            for j = 1:numel(plottedIdx)
                i = plottedIdx(j);
                yr = uniqueYears(i);
                idx = (yearVec == yr);
                [sdoy, order] = sort(doy(idx));
                tmpSyn = meanSynData(idx);
                tmpRef = meanRefData(idx);
                tmpSyn = tmpSyn(order);
                tmpRef = tmpRef(order);
                % Smooth lines
                windowSize = 15;
                tmpSyn = smoothdata(tmpSyn, 'movmean', windowSize);
                tmpRef = smoothdata(tmpRef, 'movmean', windowSize);
%                 w = 15;
%                 b = ones(1,w)/w;
%                 % apply zero-phase moving average (filtfilt)
%                 tmpSyn = filtfilt(b,1,double(tmpSyn));
%                 tmpRef = filtfilt(b,1,double(tmpRef));
                % Plot synthetic and reference
                h(j) = plot(sdoy, tmpSyn, '-', 'Color', cmap(j,:), 'LineWidth', 1.5);
                plot(sdoy, tmpRef, '--', 'Color', cmap(j,:), 'LineWidth', 1.2);
            end
            xlabel('Day of Year');
            ylabel(varLegend);
            xlim([0 366]);
            ylim([0 3.5]);
%             title(['Mean ' convertStringsToChars(targetVarL(k)) ' (All Years Overlaid)']);
            % KGE and stats for total period
%             r = corr(meanSynData, meanRefData, 'rows', 'complete');
%             alpha = std(meanSynData, 'omitnan') / std(meanRefData, 'omitnan');
%             beta = mean(meanSynData, 'omitnan') / mean(meanRefData, 'omitnan');
%             kgeSynRef = 1 - sqrt((r - 1)^2 + (alpha - 1)^2 + (beta - 1)^2);
            % Subtitle with KGE stats
%             str = {['KGE: ' num2str(kgeSynRef, '%.5f')], ...
%                 ['r: ' num2str(r, '%.5f') ', \alpha: ' num2str(alpha, '%.5f') ', \beta: ' num2str(beta, '%.5f')]} ;
%             subtitle(str);
            % Legend and aesthetics — use the synthetic handles so legend color matches both lines
            legEntries = arrayfun(@(y) num2str(y), uniqueYears(plottedIdx), 'UniformOutput', false);
            lgd = legend(h, legEntries, 'Location', 'southoutside', 'Orientation', 'Horizontal');
%             title(lgd, 'Year');
            grid on; box off;
            set(gcf, 'color', 'white');
            % Save total plot
            saveas(gcf, fullfile(outDir, ...
                ['correlation_' convertStringsToChars(targetVar(k)) '_DOYoverlay_' char(totalStart) '_' char(totalEnd) '.png']));
            close(gcf);

            % --------------------------------------------------------------------

            totalStart = startQdate(1);
            totalEnd   = endQdate(end);

            totalIndices = dates >= totalStart & dates <= totalEnd;
            totalData = validationMetric.(targetVarL(k))(totalIndices, :);

            rmseData = totalData(:,2);   % <-- RMSE time series
            selDates = dates(totalIndices);
            [yearVec, ~, ~] = ymd(selDates);
            doy = day(selDates, 'dayofyear');

            uniqueYears = unique(yearVec);
            numYears = numel(uniqueYears);
            figure('Units','inches','Position',[1 1 8 4])
            hold on;

            plottedIdx = 1:2:numYears;                 % same subsampling as before
            h = gobjects(length(plottedIdx),1);
            cmap = turbo(length(plottedIdx));

            for j = 1:numel(plottedIdx)
                i = plottedIdx(j);
                yr = uniqueYears(i);

                idx = (yearVec == yr);
                [sdoy, order] = sort(doy(idx));

                tmpRMSE = rmseData(idx);
                tmpRMSE = tmpRMSE(order);

                % Smooth (same window as before)
                windowSize = 15;
                tmpRMSE = smoothdata(tmpRMSE, 'movmean', windowSize);

                % Plot
                h(j) = plot(sdoy, tmpRMSE, '-', ...
                    'Color', cmap(j,:), 'LineWidth', 1.5);
            end
            xlabel('Day of Year');
            ylabel('RMSE [mm/day]');
            xlim([0 366]);

            % Adjust if needed
            % ylim([0 max(rmseData,[],'omitnan')*1.1]);

            legEntries = arrayfun(@(y) num2str(y), ...
                uniqueYears(plottedIdx), 'UniformOutput', false);

            % lgd = legend(h, legEntries, ...
            %     'Location', 'southoutside', ...
            %     'Orientation', 'Horizontal');

            grid on; box off;
            set(gcf, 'Color', 'white');

            saveas(gcf, fullfile(outDir, ...
                ['RMSE_' convertStringsToChars(targetVar(k)) ...
                '_DOYoverlay_' char(totalStart) '_' char(totalEnd) '.png']));

            close(gcf);

            % --------------------------------------------------------------------

            refData = refValidation.(targetVarL(k));
            synData = synImages.(targetVarL(k));
            % Handle nan values
            if ~isnan(nanValue)
                refData(refData == nanValue) = nan;
                synData(isnan(refData)) = nan;  % Mask synthetic data where reference is NaN
            end
            % === Plot variance for each individual validation period ===
            if numPeriods > 1
                for period = 1:numPeriods
                    % Time filtering
                    periodIndices = dates >= startQdate(period) & dates <= endQdate(period);
                    % Compute spatial variance over time for current period
                    varRefData = squeeze(var(refData(:,:,periodIndices), 0, [1 2], 'omitnan'));
                    varSynData = squeeze(var(synData(:,:,periodIndices), 0, [1 2], 'omitnan'));
                    % Plot
                    figure();
                    plot(dates(periodIndices), varRefData, 'r-', ...
                        dates(periodIndices), varSynData, 'k-');
                    legend('Reference', 'Synthetic', 'Location', 'northeast');
                    xlabel('Date');
                    ylabel([varLegend '^2']);
                    title([convertStringsToChars(targetVar(k)) ' variance']);
                    % KGE and stats
                    r = corr(varSynData, varRefData, 'rows', 'complete');
                    alpha = std(varSynData, 'omitnan') / std(varRefData, 'omitnan');
                    beta = mean(varSynData, 'omitnan') / mean(varRefData, 'omitnan');
                    kgeSynRef = 1 - sqrt((r - 1)^2 + (alpha - 1)^2 + (beta - 1)^2);
                    str = {['KGE: ' num2str(kgeSynRef, '%.5f')], ...
                        ['r: ' num2str(r, '%.5f') ', \alpha: ' num2str(alpha, '%.5f') ', \beta: ' num2str(beta, '%.5f')]};
                    subtitle(str);
                    grid on; box off;
                    set(gcf, 'color', 'white');
                    % Save figure
                    saveas(gcf, fullfile(outDir, ...
                        ['variance_' convertStringsToChars(targetVar(k)) '_' char(startQdate(period)) '_' char(endQdate(period)) '.png']));
                    close(gcf);
                end
            end
            % === Plot variance for the total period (from first start to last end) ===
            totalStart = startQdate(1);
            totalEnd = endQdate(end);
            totalIndices = dates >= totalStart & dates <= totalEnd;
            % Compute spatial variance over time for total period
            varRefData = squeeze(var(refData(:,:,totalIndices), 0, [1 2], 'omitnan'));
            varSynData = squeeze(var(synData(:,:,totalIndices), 0, [1 2], 'omitnan'));
            figure();
            plot(dates(totalIndices), varRefData, 'r-', ...
                dates(totalIndices), varSynData, 'k-');
            legend('Reference', 'Synthetic', 'Location', 'northeast');
            xlabel('Date');
            ylabel([varLegend '^2']);
            title([convertStringsToChars(targetVar(k)) ' variance']);
            % KGE and stats for total period
            r = corr(varSynData, varRefData, 'rows', 'complete');
            alpha = std(varSynData, 'omitnan') / std(varRefData, 'omitnan');
            beta = mean(varSynData, 'omitnan') / mean(varRefData, 'omitnan');
            kgeSynRef = 1 - sqrt((r - 1)^2 + (alpha - 1)^2 + (beta - 1)^2);
            str = {['KGE: ' num2str(kgeSynRef, '%.5f')], ...
                ['r: ' num2str(r, '%.5f') ', \alpha: ' num2str(alpha, '%.5f') ', \beta: ' num2str(beta, '%.5f')]};
            subtitle(str);
            grid on; box off;
            set(gcf, 'color', 'white');
            % Save total plot
            saveas(gcf, fullfile(outDir, ...
                ['variance_' convertStringsToChars(targetVar(k)) '_total_' char(totalStart) '_' char(totalEnd) '.png']));
            close(gcf);

            % -------------------------------------------------------------------------

            refData = refValidation.(targetVarL(k));
            synData = synImages.(targetVarL(k));
            % Handle nan values
            if ~isnan(nanValue)
                refData(refData == nanValue) = nan;
                synData(isnan(refData)) = nan;
            end
            % === Loop over each validation period ===
            if numPeriods > 1
                for period = 1:numPeriods
                    % Time slice
                    periodIndices = dates >= startQdate(period) & dates <= endQdate(period);
                    refSlice = refData(:,:,periodIndices);
                    synSlice = synData(:,:,periodIndices);
                    % Initialize error matrices
                    absDayErr = abs(synSlice - refSlice);
                    dayErr    = synSlice - refSlice;
                    relErr    = absDayErr ./ refSlice;
                    % Clean up
                    relErr(relErr == Inf | relErr == -Inf) = nan;   % <--------------------------------------------------------------
                    % Temporal average of error metrics
                    meanError  = mean(absDayErr, 3, 'omitnan');
                    meanBias   = mean(dayErr, 3, 'omitnan');
                    meanRelErr = mean(relErr, 3, 'omitnan');

                    % MAE
                    figure;
                    figMean = imagesc(meanError);
                    set(figMean, 'AlphaData', ~isnan(meanError));
                    colormap(gca, turbo(256));
                    title('Mean absolute error');
                    subtitle(['Mean MAE: ' num2str(mean(mean(meanError, 'omitnan'), 'omitnan'), '%1.5f')]);
                    set(gcf, 'color', 'white');
                    hcb = colorbar;
                    set(get(hcb, 'label'), 'string', 'Mean absolute error [mm/day]', 'Rotation', 90);
                    axis equal off;
                    saveas(gcf, fullfile(outDir, ...
                        ['mae_' convertStringsToChars(targetVar(k)) '_' char(startQdate(period)) '_' char(endQdate(period)) '.png']));
                    close(gcf);

                    % Bias
                    figure;
                    figMean = imagesc(meanBias);
                    set(figMean, 'AlphaData', ~isnan(meanBias));
                    colormap(gca, turbo(256));
                    title('Bias');
                    subtitle(['Mean bias: ' num2str(mean(mean(meanBias, 'omitnan'), 'omitnan'), '%1.5f')]);
                    set(gcf, 'color', 'white');
                    hcb = colorbar;
                    set(get(hcb, 'label'), 'string', 'Bias [mm/day]', 'Rotation', 90);
                    axis equal off;
                    saveas(gcf, fullfile(outDir, ...
                        ['bias_' convertStringsToChars(targetVar(k)) '_' char(startQdate(period)) '_' char(endQdate(period)) '.png']));
                    close(gcf);

                    % MRE
                    figure;
                    figMean = imagesc(meanRelErr);
                    set(figMean, 'AlphaData', ~isnan(meanRelErr));
                    colormap(gca, turbo(256));
                    title('Mean relative error');
                    subtitle(['Mean: ' num2str(mean(mean(meanRelErr, 'omitnan'), 'omitnan'), '%1.5f')]);
                    set(gcf, 'color', 'white');
                    hcb = colorbar;
                    set(get(hcb, 'label'), 'string', 'Mean relative error', 'Rotation', 90);
                    axis equal off;
                    saveas(gcf, fullfile(outDir, ...
                        ['mre_' convertStringsToChars(targetVar(k)) '_' char(startQdate(period)) '_' char(endQdate(period)) '.png']));
                    close(gcf);
                end
            end
            % === Now process the total period ===
            totalStart = startQdate(1);
            totalEnd   = endQdate(end);
            totalIndices = dates >= totalStart & dates <= totalEnd;
            refSlice = refData(:,:,totalIndices);
            synSlice = synData(:,:,totalIndices);
            absDayErr = abs(synSlice - refSlice);
            dayErr    = synSlice - refSlice;
            relErr    = absDayErr ./ refSlice;
            relErr(relErr == Inf | relErr == -Inf) = nan; % <-----------------------------------------------------------------------------------
            meanError  = mean(absDayErr, 3, 'omitnan');
            meanBias   = mean(dayErr, 3, 'omitnan');
            meanRelErr = mean(relErr, 3, 'omitnan');

            % MAE (Total)
            figure('Units','inches','Position',[1 1 3 3])
            figMean = imagesc(meanError);
            set(figMean, 'AlphaData', ~isnan(meanError));
            colormap(gca, turbo(256));
            % title('Mean absolute error');
            % subtitle(['Mean MAE: ' num2str(mean(mean(meanError, 'omitnan'), 'omitnan'), '%1.5f')]);
            set(gcf, 'color', 'white');
            hcb = colorbar;
            hcb.FontSize = 12;
            hcb.Label.String = 'MAE [mm/day]';
            hcb.Label.FontSize = 12;
            hcb.Label.Rotation = 90;
            axis equal off;
            saveas(gcf, fullfile(outDir, ...
                ['mae_' convertStringsToChars(targetVar(k)) '_total_' char(totalStart) '_' char(totalEnd) '.png']));
            close(gcf);

            % Bias (Total)
            figure('Units','inches','Position',[1 1 3 3])
            figMean = imagesc(meanBias);
            set(figMean, 'AlphaData', ~isnan(meanBias));
            colormap(gca, turbo(256));
            % title('Bias');
            % subtitle(['Mean bias: ' num2str(mean(mean(meanBias, 'omitnan'), 'omitnan'), '%1.5f')]);
            set(gcf, 'color', 'white');
            hcb = colorbar;
            hcb.FontSize = 12;
            hcb.Label.String = 'Bias [mm/day]';
            hcb.Label.FontSize = 12;
            hcb.Label.Rotation = 90;
            hcb.Label.Position(1) = hcb.Label.Position(1) - 0.8;
            axis equal off;
            saveas(gcf, fullfile(outDir, ...
                ['bias_' convertStringsToChars(targetVar(k)) '_total_' char(totalStart) '_' char(totalEnd) '.png']));
            close(gcf);

            % MRE (Total)
            figure('Units','inches','Position',[1 1 3 3])
            figMean = imagesc(meanRelErr);
            set(figMean, 'AlphaData', ~isnan(meanRelErr));
            colormap(gca, turbo(256));
            % title('Mean relative error');
            % subtitle(['Mean: ' num2str(mean(mean(meanRelErr, 'omitnan'), 'omitnan'), '%1.5f')]);
            set(gcf, 'color', 'white');
            hcb = colorbar;
            hcb.FontSize = 12;
            hcb.Label.String = 'MRE [-]';
            hcb.Label.FontSize = 12;
            hcb.Label.Rotation = 90;
            axis equal off;
            saveas(gcf, fullfile(outDir, ...
                ['mre_' convertStringsToChars(targetVar(k)) '_total_' char(totalStart) '_' char(totalEnd) '.png']));
            close(gcf);

            % -------------------------------------------------------------------------

            targetVarL = lower(targetVar(k));
            for cvVar = 1:numel(climateVars)
                % Extract covariate data
                covName = lower(climateVars(cvVar));
                cov = table2array(climateData(:,covName));
                cv = cat(3, cov{:});
                datesCov = datetime(table2array(climateData(:,"date")),'ConvertFrom','yyyyMMdd','Format','dd-MM-yyyy');

                datasets = struct('name', {'Synthetic', 'Real'}, 'data', {synImages, refValidation});

                for period = 1:numel(startQdate)
                    densities = cell(1, numel(datasets));
                    histograms = cell(1, numel(datasets));
                    binCentersX = [];
                    binCentersY = [];

                    % Loop through syn and ref
                    for d = 1:numel(datasets)
                        tg = datasets(d).data.(targetVarL);
                        datesTar = datetime(datasets(d).data.date, 'ConvertFrom','yyyyMMdd','Format','dd-MM-yyyy');

                        idxTar = datesTar >= startQdate(period) & datesTar <= endQdate(period);
                        idxCov = datesCov >= startQdate(period) & datesCov <= endQdate(period);

                        x = tg(:,:,idxTar);
                        y = cv(:,:,idxCov);
                        x = x(:);
                        y = y(:);

                        validIdx = isfinite(x) & isfinite(y);
                        x = x(validIdx);
                        y = y(validIdx);

                        % Define bins
                        nbins = 100;
                        edgesX = linspace(min(x), max(x), nbins);
                        edgesY = linspace(min(y), max(y), nbins);
                        binCentersX = edgesX(1:end-1) + diff(edgesX)/2;
                        binCentersY = edgesY(1:end-1) + diff(edgesY)/2;

                        binX = discretize(x, edgesX);
                        binY = discretize(y, edgesY);

                        validBins = ~isnan(binX) & ~isnan(binY);
                        xPlot = x(validBins);
                        yPlot = y(validBins);
                        binX = binX(validBins);
                        binY = binY(validBins);

                        binLinear = sub2ind([nbins, nbins], binX, binY);
                        binCounts = accumarray(binLinear, 1);
                        pointDensity = binCounts(binLinear);

                        densities{d} = pointDensity;

                        % Compute 2D histogram
                        hist2D = accumarray([binX, binY], 1, [nbins, nbins]);
                        histograms{d} = hist2D;
                    end

                    %% Match color scale
                    allDensities = vertcat(densities{:});
                    clim = [min(allDensities), max(allDensities)];

                    % Plot syn and ref
                    for d = 1:numel(datasets)
                        tg = datasets(d).data.(targetVarL);
                        datesTar = datetime(datasets(d).data.date, 'ConvertFrom','yyyyMMdd','Format','dd-MM-yyyy');

                        idxTar = datesTar >= startQdate(period) & datesTar <= endQdate(period);
                        idxCov = datesCov >= startQdate(period) & datesCov <= endQdate(period);

                        x = tg(:,:,idxTar);
                        y = cv(:,:,idxCov);
                        x = x(:);
                        y = y(:);
                        validIdx = isfinite(x) & isfinite(y);
                        x = x(validIdx);
                        y = y(validIdx);

                        edgesX = linspace(min(x), max(x), nbins);
                        edgesY = linspace(min(y), max(y), nbins);
                        binX = discretize(x, edgesX);
                        binY = discretize(y, edgesY);
                        validBins = ~isnan(binX) & ~isnan(binY);
                        xPlot = x(validBins);
                        yPlot = y(validBins);
                        binX = binX(validBins);
                        binY = binY(validBins);

                        binLinear = sub2ind([nbins, nbins], binX, binY);
                        binCounts = accumarray(binLinear, 1);
                        pointDensity = binCounts(binLinear);

                        figure('Color','white','Position',[100 100 700 600]);
                        scatter(xPlot, yPlot, 5, pointDensity, 'filled', 'MarkerFaceAlpha', 0.5);
                        colormap(turbo);
                        caxis(clim);  % same colorbar limits
                        colorbar;
                        xlabel('Et');
                        ylabel(covName);
                        title(sprintf('%s – %s to %s', datasets(d).name, datestr(startQdate(period)), datestr(endQdate(period))));
                        grid on;
                        axis tight;

                        fname = sprintf('%s_densityScatter_%s_%s_%s_%s.png', ...
                            datasets(d).name, targetVarL, covName, ...
                            char(startQdate(period)), char(endQdate(period)));
                        saveas(gcf, fullfile(outDir, fname));
                        close(gcf);
                    end

                    %% Compute histogram-based metrics
                    % Normalize histograms
                    h1 = histograms{1} / sum(histograms{1}(:));
                    h2 = histograms{2} / sum(histograms{2}(:));

                    % RMSE
                    rmse = sqrt(mean((h1(:) - h2(:)).^2));

                    % Hellinger Distance
                    hellinger = sqrt(0.5 * sum((sqrt(h1(:)) - sqrt(h2(:))).^2));

                    fprintf('Period %d [%s to %s]:\n', period, ...
                        datestr(startQdate(period)), datestr(endQdate(period)));
                    fprintf('  RMSE = %.4f\n', rmse);
                    fprintf('  Hellinger = %.4f\n\n', hellinger);
                end
            end

            % -------------------------------------------------------------------------

            if outputTime == 1 && pixelWise == false
                cellData = synImages.(strcat(targetVarL(k), "_Distances"));
                minValues    = cellfun(@min, cellData);
                maxValues    = cellfun(@max, cellData);
                medianValues = cellfun(@median, cellData);
                numK = numel(cellData{1});

                % === Loop over each validation period ===
                if numPeriods > 1
                    for period = 1:numPeriods
                        % Period selection
                        idx = dates >= startQdate(period) & dates <= endQdate(period);
                        periodDates      = dates(idx);
                        periodMin        = minValues(idx);
                        periodMax        = maxValues(idx);
                        periodMedian     = medianValues(idx);
                        periodCellData   = cellData(idx);
                        % === Plot for the current period ===
                        figure; hold on;
                        fill([periodDates; flipud(periodDates)], [periodMin; flipud(periodMax)], ...
                            'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                        for i = 1:numel(periodCellData)
                            scatter(repmat(periodDates(i), numK, 1), periodCellData{i}, ...
                                20, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
                        end
                        plot(periodDates, periodMedian, 'r--', 'LineWidth', 1.5);
                        xlabel('Date');
                        switch metricKNN
                            case 1, ylabel('RMSE');
                            case 2, ylabel('MAE');
                            case 3, ylabel('1 - bSPEM');
                            case 4, ylabel('Hellinger Distance');
                            case 5, ylabel('(1 - bSPEM) + Hellinger');
                            case 6, ylabel('SPAEF');
                        end
                        title('Daily distance of the k candidates');
                        grid on; hold off;
                        set(gcf, 'color', 'white');
                        % Save figure with period in filename
                        saveas(gcf, fullfile(outDir, ...
                            ['distance_' convertStringsToChars(targetVar(k)) '_' ...
                            char(startQdate(period)) '_' char(endQdate(period)) '.png']));
                        close(gcf);
                    end
                end
                % === Plot for the full validation period ===
                totalStart = startQdate(1);
                totalEnd   = endQdate(end);
                idx = dates >= totalStart & dates <= totalEnd;
                totalDates      = dates(idx);
                totalMin        = minValues(idx);
                totalMax        = maxValues(idx);
                totalMedian     = medianValues(idx);
                totalCellData   = cellData(idx);
                figure; hold on;
                fill([totalDates; flipud(totalDates)], [totalMin; flipud(totalMax)], ...
                    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                for i = 1:numel(totalCellData)
                    scatter(repmat(totalDates(i), numK, 1), totalCellData{i}, ...
                        20, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
                end
                plot(totalDates, totalMedian, 'r--', 'LineWidth', 1.5);
                xlabel('Date');
                switch metricKNN
                    case 1, ylabel('RMSE');
                    case 2, ylabel('MAE');
                    case 3, ylabel('1 - bSPEM');
                    case 4, ylabel('Hellinger Distance');
                    case 5, ylabel('(1 - bSPEM) + Hellinger');
                    case 6, ylabel('SPAEF');
                end
                title(['Daily distance of the k candidates (total period)']);
                grid on; hold off;
                set(gcf, 'color', 'white');
                % Save figure for total period
                saveas(gcf, fullfile(outDir, ...
                    ['distance_' convertStringsToChars(targetVar(k)) '_total_' ...
                    char(totalStart) '_' char(totalEnd) '.png']));
                close(gcf);
            end

            % -------------------------------------------------------------------------

            if createGIF == true
                % Set the output GIF file name
                if pixelWise == false
                    gifVal = fullfile(outDir,['\validation_' convertStringsToChars(targetVar(k)) '.gif']);

                    refDates = refValidation.date;
                    synDates = synImages.date;
                    dates    = datetime(synDates,'ConvertFrom','yyyyMMdd','format','dd/MM/yyyy');
                    bdName   = [convertStringsToChars(targetVar(k)) '_BestDistance'];
                    analogs  = sortedDates(:,2);
                    currentDOY  = nan(numel(analogs{1}),numel(dates));
                    meanDOY     = nan(numel(dates),1);
                    diffDOY     = nan(numel(analogs{1}),1);
                    diffBest    = nan(numel(dates),1);
                    %bestAnalog  = nan(numel(dates),1);
                    %bestDist    = synImages.(bdName);
                    %currentBest = nan(size(synDates));
                    q1DOY     = [];
                    q3DOY     = [];
                    minOut      = [];
                    maxOut      = [];
                    minCurDOY   = nan(numel(dates),1);
                    maxCurDOY   = nan(numel(dates),1);

                    % Create an empty figure
                    figure();
                    % Loop over each file in the synthetic directory and find the corresponding
                    % file in the reference directory
                    for i = 1:size(synData,3)
                        % Find the corresponding file in the reference directory with the same name
                        referenceIndex = find(refDates == synDates(i));

                        analogDOY  = day(datetime(sort(analogs{i},'descend'),'ConvertFrom','yyyyMMdd'),'dayofyear');
                        bestAnalog = day(datetime(analogs{i}(1),'ConvertFrom','yyyyMMdd'),'dayofyear');
                        refDOY     = day(datetime(refDates(i),'ConvertFrom','yyyyMMdd'),'dayofyear');
                        % Check for circular transition
                        % Compute the minimum and maximum range boundaries
                        minRange = refDOY - daysRange;
                        maxRange = refDOY + daysRange;
                        % Handle circular transition for minRangeQ
                        if minRange <= 0
                            minRange = 365 + minRange;
                        end
                        % Handle circular transition for maxRangeQ
                        if maxRange > 365
                            maxRange = maxRange - 365;
                        end
                        % Construct the rangeTot array
                        if minRange < maxRange
                            rangeTot = minRange:maxRange;
                        else
                            rangeQmin = minRange:365;
                            rangeQmax = 1:maxRange;
                            rangeTot    = [rangeQmin rangeQmax];
                        end
                        if refDOY == 366
                            refDOY = 1;
                        end
                        rangeTot = rangeTot';
                        % Find the index of refDOY and analogDOY(i) in rangeTot
                        indexRef = find(rangeTot == refDOY);
                        for di = 1:numel(analogDOY)
                            if analogDOY(di) == 366
                                analogDOY(di) = 1;
                            end
                            indexAnalog   = find(rangeTot == analogDOY(di));
                            diffDOY(di)   = indexRef - indexAnalog;
                        end
                        if bestAnalog == 366
                            bestAnalog = 1;
                        end
                        idxBestAnalog = find(rangeTot == bestAnalog);
                        diffBest(i)   = indexRef - idxBestAnalog;

                        % If a matching file is found, display the two images side by side
                        if ~isempty(referenceIndex)
                            % Load the two images
                            synthetic = synData(:,:,i);
                            %synthetic(synthetic==min(min(synthetic))) = NaN;
                            %synthetic(synthetic==-999) = NaN;
                            reference = refData(:,:,referenceIndex);
                            %reference(reference==min(min(reference))) = NaN;
                            %reference(reference==-999) = NaN;

                            sgtitle(targetVar(k))

                            % Create a figure with three subplots
                            subplot(3,3,[1,4]);
                            img1 = imshow(synthetic);
                            colormap(gca, turbo(256));
                            set(img1, 'AlphaData', ~isnan(synthetic))
                            %caxis([0 maxColor])
                            caxis(varRange)
                            axis equal
                            title('Synthetic');
                            %colorbar(gca,'southoutside')

                            subplot(3,3,[2,5]);
                            img2 = imshow(reference);
                            colormap(gca, turbo(256));
                            set(img2, 'AlphaData', ~isnan(synthetic))
                            %caxis([0 maxColor])
                            caxis(varRange)
                            axis equal
                            title('Reference');
                            h = colorbar(gca,'southoutside');

                            % ERROR MAP
                            error  = synthetic - reference;
                            subplot(3,3,[3,6])
                            errMap = imshow(error);
                            set(errMap, 'AlphaData', ~isnan(synthetic))
                            colormap(gca, coolwarm(256));
                            caxis(errRange)
                            title('Error');
                            axis equal
                            h_err = colorbar(gca,'southoutside');

                            % Add a colorbar to the reference image subplot
                            %h = colorbar('southoutside');
                            set(h, 'Position', [0.13 0.4 0.5 0.03]);
                            set(get(h,'label'),'string',varLegend);
                            set(h_err, 'Position', [0.7 0.4 0.205 0.03])
                            set(get(h_err,'label'),'string',varLegend);

                            % DOY difference
                            subplot(3,3,[7,8,9])
                            currentDOY(:,i) = diffDOY;
                            meanDOY(i)      = mean(currentDOY(:,i));
                            %if numel(dates)<60
                            %    boxchart(currentDOY);%,synDates(i));
                            %else
                            if nbImages > 1
                                inBetweenRegionX = [(1:i), fliplr(1:i)];
                                %inBetweenRegionX = [1:numel(dates), fliplr(1:numel(dates))];
                                quartileDOY = quantile(currentDOY(:,i),4);
                                q1DOY = [q1DOY quartileDOY(1)];
                                q3DOY = [q3DOY quartileDOY(3)];
                                inBetweenRegionY = [q3DOY, fliplr(q1DOY)];
                                itqDOY(i) = iqr(currentDOY(:,i));
                                minOut = [minOut, min(currentDOY(currentDOY(:,i) > (q1DOY(i) - itqDOY(i)),i))];
                                maxOut = [maxOut, max(currentDOY(currentDOY(:,i) < (q3DOY(i) + itqDOY(i)),i))];
                                inBetweenRegionY2 = [minOut, fliplr(maxOut)];
                                minCurDOY(i) = min(currentDOY(:,i));
                                maxCurDOY(i) = max(currentDOY(:,i));
                                patch(inBetweenRegionX, inBetweenRegionY, 'k', 'LineStyle', 'none', 'FaceAlpha', 0.15)
                                hold on
                                patch(inBetweenRegionX, inBetweenRegionY2, 'k', 'LineStyle', 'none', 'FaceAlpha', 0.1)
                                hold on
                                plot(1:i,minCurDOY(1:i),'LineStyle','--','Color','k')
                                hold on
                                plot(1:i,maxCurDOY(1:i),'LineStyle','--','Color','k')
                            end
                            %end
                            hold on
                            if nbImages > 1
                                plot(1:i,meanDOY(1:i),"Color",'red')
                            end
                            scatter(1:i,diffBest(1:i),"red")
                            %plot(dates(1:i),meanDOY(1:i),"Color",'red')
                            hold off
                            xlim([1 numel(dates)])
                            %xlim([min(dates) max(dates)])
                            ylim([-daysRange-10 daysRange+10])
                            title(['Mean DOY difference: ' num2str(mean(diffDOY),'%2.0f')]);
                            ylabel('DOY difference')
                            xlabel('Date')
                            grid on
                            ax = gca();
                            %ax.XTick = categorical(1:0.5:numel(dates));
                            %ax.XTickLabels = char(dates);

                            %                 % Best candidate MAE
                            %                 subplot(3,3,[7,8,9])
                            %                 currentBest(1:i) = bestDist(1:i);
                            %                 plot(dates,currentBest);
                            %                 hold on
                            %                 plot(dates(i),bestDist(i),"Marker","o","Color",'red')
                            %                 hold off
                            %                 xlim([min(dates) max(dates)])
                            %                 ylim([0.5 2])
                            %                 title(['Best candidate MAE: ' num2str(bestDist(i),'%1.5f')]);
                            %                 ylabel('MAE')
                            %                 xlabel('Date')
                            %                 grid on

                            % Set the title of the figure to the name of the images
                            if metricV == 1
                                sgtitle({['{\bf\fontsize{14}' char(dates(i)) '}'], ...
                                    ['{\fontsize{13}' 'MAE: ' num2str(validationMetric.(targetVarL(k))(i,2),'%1.5f') '}']})
                            elseif metricV == 2
                                sgtitle({['{\bf\fontsize{14}' char(dates(i)) '}'], ...
                                    ['{\fontsize{13}' 'RMSE: ' num2str(validationMetric.(targetVarL(k))(i,2),'%1.5f') '}']})
                            elseif metricV == 3
                                sgtitle({['{\bf\fontsize{14}' char(dates(i)) '}'], ...
                                    ['{\fontsize{13}' 'SPEM: ' num2str(validationMetric.(targetVarL(k))(i,2),'%1.5f') '}']})
                            elseif metricV == 4
                                sgtitle({['{\bf\fontsize{14}' char(dates(i)) '}'], ...
                                    ['{\fontsize{13}' 'SPAEF: ' num2str(validationMetric.(targetVarL(k))(i,2),'%1.5f') '}']})
                            elseif metricV == 5
                                sgtitle({['{\bf\fontsize{14}' char(dates(i)) '}'], ...
                                    ['{\fontsize{13}' 'KGE: ' num2str(validationMetric.(targetVarL(k))(i,2),'%1.5f') '}']})
                            end

                            % Save the current frame as a GIF
                            set(gcf, 'color', 'white');
                            frame = getframe(gcf);
                            im = frame2im(frame);
                            [imind, cm] = rgb2ind(im, 256);
                            if i == 1
                                imwrite(imind, cm, gifVal, 'gif', 'Loopcount', size(synData,3), 'DelayTime', 0.1);
                            else
                                imwrite(imind, cm, gifVal, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
                            end
                        end
                    end
                else

                    % -------------------------------------------------------------------------

                    % Set the output GIF file name
                    gifVal = fullfile(outDir,['\validation_' convertStringsToChars(targetVar(k)) '.gif']);

                    refDates = refValidation.date;
                    refData  = refValidation.(targetVarL(k));
                    synDates = synImages.date;
                    synData  = synImages.(targetVarL(k));
                    dates    = datetime(synDates,'ConvertFrom','yyyyMMdd','format','dd/MM/yyyy');

                    % Create an empty figure
                    figure();
                    % Loop over each file in the synthetic directory and find the corresponding
                    % file in the reference directory
                    for i = 1:size(synData,3)
                        % Find the corresponding file in the reference directory with the same name
                        referenceIndex = find(refDates == synDates(i));

                        % If a matching file is found, display the two images side by side
                        if ~isempty(referenceIndex)
                            % Load the two images
                            synthetic = synData(:,:,i);
                            %synthetic(synthetic==min(min(synthetic))) = NaN;
                            %synthetic(synthetic==-999) = NaN;
                            reference = refData(:,:,referenceIndex);
                            %reference(reference==min(min(reference))) = NaN;
                            %reference(reference==-999) = NaN;

                            sgtitle(targetVar(k))

                            % Create a figure with three subplots
                            subplot(1,3,1);
                            img1 = imshow(synthetic);
                            colormap(gca, turbo(256));
                            set(img1, 'AlphaData', ~isnan(synthetic))
                            %caxis([0 maxColor])
                            caxis(varRange)
                            axis equal
                            title('Synthetic');
                            %colorbar(gca,'southoutside')

                            subplot(1,3,2);
                            img2 = imshow(reference);
                            colormap(gca, turbo(256));
                            set(img2, 'AlphaData', ~isnan(synthetic))
                            %caxis([0 maxColor])
                            caxis(varRange)
                            axis equal
                            title('Reference');
                            h = colorbar(gca,'southoutside');

                            % ERROR MAP
                            error  = synthetic - reference;
                            subplot(1,3,3)
                            errMap = imshow(error);
                            set(errMap, 'AlphaData', ~isnan(synthetic))
                            colormap(gca, coolwarm(256));
                            caxis(errRange)
                            title('Error');
                            axis equal
                            h_err = colorbar(gca,'southoutside');

                            % Add a colorbar to the reference image subplot
                            %h = colorbar('southoutside');
                            set(h, 'Position', [0.13 0.2 0.5 0.03]);
                            set(get(h,'label'),'string',varLegend);
                            set(h_err, 'Position', [0.7 0.2 0.205 0.03])
                            set(get(h_err,'label'),'string',varLegend);

                            % Set the title of the figure to the name of the images
                            if metricV == 1
                                sgtitle({['{\bf\fontsize{14}' char(dates(i)) '}'], ...
                                    ['{\fontsize{13}' 'MAE: ' num2str(validationMetric.(targetVarL(k))(i,2),'%1.5f') '}']})
                            elseif metricV == 2
                                sgtitle({['{\bf\fontsize{14}' char(dates(i)) '}'], ...
                                    ['{\fontsize{13}' 'RMSE: ' num2str(validationMetric.(targetVarL(k))(i,2),'%1.5f') '}']})
                            elseif metricV == 3
                                sgtitle({['{\bf\fontsize{14}' char(dates(i)) '}'], ...
                                    ['{\fontsize{13}' 'SPEM: ' num2str(validationMetric.(targetVarL(k))(i,2),'%1.5f') '}']})
                            elseif metricV == 4
                                sgtitle({['{\bf\fontsize{14}' char(dates(i)) '}'], ...
                                    ['{\fontsize{13}' 'SPAEF: ' num2str(validationMetric.(targetVarL(k))(i,2),'%1.5f') '}']})
                            elseif metricV == 5
                                sgtitle({['{\bf\fontsize{14}' char(dates(i)) '}'], ...
                                    ['{\fontsize{13}' 'KGE: ' num2str(validationMetric.(targetVarL(k))(i,2),'%1.5f') '}']})
                            end

                            % Save the current frame as a GIF
                            set(gcf, 'color', 'white');
                            frame = getframe(gcf);
                            im = frame2im(frame);
                            [imind, cm] = rgb2ind(im, 256);
                            if i == 1
                                imwrite(imind, cm, gifVal, 'gif', 'Loopcount', size(synData,3), 'DelayTime', 0.1);
                            else
                                imwrite(imind, cm, gifVal, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
                            end
                        end
                    end
                end
            end
        end
    end
end

end
