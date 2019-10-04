function [Energy_Resolution, res_fit, fitresult] = Reconstructed_Energy(output,hist_channel_max,Filt,filtering,depict,en_resol)
% Reconstructed Energy Histogram
% permette di mostrare l'istogramma di energia (con i filtri desiderati
% applicati) espresso sia in canali ADC che calibrato in energia e restituisce il valore di risoluzione energetica ottenuto

constant_FWHM = 2.355;
n_bins = 5000;
space = hist_channel_max/n_bins;
x_ch = [0:space:hist_channel_max];

% FILTERING ON RECONSTRUCTED EVENTS: Events are filtered according to reconstruction error and spatial limits of recontructed positions
 if filtering == 1
        % Filtro gli eventi sulla base dell'errore di ricostruzione e della
        % regione spaziale del cristallo di interesse
        F3=(output.error>=Filt.error_min).*(output.error<=Filt.error_max);
        filt_x = (output.x_rec>Filt.x_rec_min).*(output.x_rec<Filt.x_rec_max);
        filt_y = (output.y_rec>Filt.y_rec_min).*(output.y_rec<Filt.y_rec_max);
        F4 = (filt_x.*filt_y);

        % Search for the indices of the events that match al the set filters
        Filt.kill_ML = find(F3.*F4);
        
        %Events that pass ALL the set filters are selected
        Energy = output.energy(Filt.kill_ML);
        
 else 
         Energy = output.energy;
 end
 
 % SELECTION OF ENERGY SPECTRUM FOR FITTING: The energy spectrum[ADC channels] of all the events passing the set
 % filters is plotted and the user has to manually select the histogram
 % portion for gaussian fitting
 if en_resol == 1
    figure
    hist(Energy, 5000);
    title('SELECT THE SPECTRUM PORTION FOR GAUSSIAN FITTING','Fontname','Arial','FontWeight','bold')
    [Channels] = ginput(2);
    Energy = Energy(find((Energy>Channels(1,1)).*(Energy< Channels(2,1))));
    
 end
 [Counts, Centers]=hist(Energy,x_ch);
 
%The selected histogram, the gaussian fitting will be applied to, is plotted: on the x-axis the energy amplitude is expressed
%in ADC-channels

% Reconstructed energy spectrum [ADC channels]
figure
plot(Centers, Counts ,'LineWidth',2.5)
xlim([prctile(Energy,0.1) prctile(Energy,99.9)])
ylim([0 1.1*max(Counts)])
hold on

set(gca,'Fontsize',14,'Fontname','Arial','FontWeight','bold')
xlabel('Reconstructed Energy [ADC channels]','Fontname','Arial','FontWeight','bold')
ylabel('Counts [a.u]','Fontname','Arial','FontWeight','bold')
title('Reconstructed Signal Spectrum (Statistical Method)','Fontname','Arial','FontWeight','bold')
grid on

if en_resol ==1 
    %ENERGY CALIBRATION: the x-axis of the spectrum is converted from ADC
    %channels into energy values
    fitresult = fit(x_ch', Counts','gauss2'); % the gaussian fitting is used to search for the peaks in the spectrum
    coeffvals = coeffvalues(fitresult)';
    coeffvals = reshape(coeffvals,size(coeffvals,1)/2,2);
    mu = sort(coeffvals(2,:));
    %Alternativa: selezione manuale dei picchi degli spettri
    % [Peaks] = ginput(2);
    % mu = sort(Peaks(:,1))';

    E = [Filt.E_calibration_min Filt.E_calibration_max];
    res_fit = fit(mu',E','poly1');
    calibration_line = coeffvalues(res_fit);
    x_en = res_fit(Centers); % x-axis calibration: energy values in ADC channels are converted into keV energy values

    if depict == 1
        % Energy calibration line plot
        x_plot = linspace(0,mu(2)*1.4, 1000);
        figure
        plot(mu,E,'o')
        hold on
        plot(x_plot ,res_fit(x_plot),'--')
        title(strcat('Energy calibration line: m = ', num2str(calibration_line(1)),'; q = ', num2str(calibration_line(2))))
        xlabel('Energy [ADC channels or a.u.]')
        ylabel('Energy [keV]')
         legend('Calibration point','Fitting')

        % Calibrated energy spectrum plot
        figure
        plot(x_en,Counts)
        title('Calibrated Energy Spectrum')
        xlabel('Energy [keV]')
        ylabel('Counts [a.u.]')
        hold on
    end


    % ENERGY RESOLUTION COMPUTATION: Compute energy resolution
    fitresult1 = fit(x_en, Counts','gauss2');
    coeffvals = coeffvalues(fitresult1)';
    coeffvals = reshape(coeffvals,size(coeffvals,1)/2,2);
    photopeak_index = find(max(coeffvals(2,:)));
    FWHM_photopeak = constant_FWHM *coeffvals(3,photopeak_index)/sqrt(2);
    Energy_Resolution = FWHM_photopeak/(coeffvals(2,photopeak_index))*100;


    if depict == 1

        % Gaussian fitting of the calibrated energy spectrum plot
        plot(x_en,fitresult1(x_en),'--')
        title('Gaussian fitting of the Calibrated Energy Spectrum')
        legend('Spectrum','Fitting')
        xlabel('Energy [keV]')
        ylabel('Counts [a.u.]')

    end

else
    Energy_Resolution = NaN;
    res_fit = NaN;
    fitresult = NaN;
end

end