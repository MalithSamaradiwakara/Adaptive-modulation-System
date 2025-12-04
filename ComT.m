clc;
clear;
close;

% Simulation Parameters
N = 1e6;
N = N - mod(N,4);
SNRdB_range = 0:2:20;
channelType = 'Rician';
ricianK = 5;
numFrames = 50;

schemes = {'BPSK','QPSK','16QAM'};
bitsPerSymbol = [1,2,4];

% Preallocate result matrices
BER_results = zeros(length(schemes), length(SNRdB_range));
SE_results = zeros(length(schemes), length(SNRdB_range));
BER_adaptive = zeros(1,length(SNRdB_range));
SE_adaptive = zeros(1,length(SNRdB_range));
modulationChoice = zeros(1,length(SNRdB_range));
targetBER = 1e-3;

% Modulation Simulation
for s = 1:length(schemes)
    scheme = schemes{s};
    k = bitsPerSymbol(s);
    fprintf("\nRunning %s...\n", scheme);

    for snrIdx = 1:length(SNRdB_range)
        SNRdB = SNRdB_range(snrIdx);
        totalErrors = 0; totalBits = 0;

        for frame = 1:numFrames
            bits = randi([0 1], N, 1);
            symbols = modulate_bits(bits, scheme);
            rxSymbols = channel_model(symbols, SNRdB, channelType, ricianK);
            rxBits = demodulate_symbols(rxSymbols, scheme);
            % Count errors
            totalErrors = totalErrors + sum(bits ~= rxBits);
            totalBits = totalBits + length(bits);
        end
        % Store BER and Spectral Efficiency
        BER_results(s, snrIdx) = totalErrors / totalBits;
        SE_results(s, snrIdx) = k * (1 - BER_results(s, snrIdx));
    end
end

% Adaptive Modulation Simulation
fprintf("\nRunning Adaptive Modulation (BER-target based)...\n");
for snrIdx = 1:length(SNRdB_range)
    SNRdB = SNRdB_range(snrIdx);
    SNRlin = 10^(SNRdB/10);

    %Calculate theoretical BER for each scheme
    BER_BPSK  = qfunc(sqrt(2*SNRlin));
    BER_QPSK  = qfunc(sqrt(2*SNRlin));
    BER_16QAM = 3/8 * erfc(sqrt(0.1*SNRlin));

    %Choose highest order modulation that meets target BER
    if BER_16QAM < targetBER
        adaptiveScheme = '16QAM'; k = 4; modulationChoice(snrIdx) = 4;
    elseif BER_QPSK < targetBER
        adaptiveScheme = 'QPSK'; k = 2; modulationChoice(snrIdx) = 2;
    else
        adaptiveScheme = 'BPSK'; k = 1; modulationChoice(snrIdx) = 1;
    end

    %Run simulation for adaptive modulation
    totalErrors = 0; totalBits = 0;
    for frame = 1:numFrames
        bits = randi([0 1], N, 1);
        symbols = modulate_bits(bits, adaptiveScheme);
        rxSymbols = channel_model(symbols, SNRdB, channelType, ricianK);
        rxBits = demodulate_symbols(rxSymbols, adaptiveScheme);

        totalErrors = totalErrors + sum(bits ~= rxBits);
        totalBits = totalBits + length(bits);
    end

    % Store BER and Spectral Efficiency
    BER_adaptive(snrIdx) = totalErrors / totalBits;
    SE_adaptive(snrIdx) = k * (1 - BER_adaptive(snrIdx));
end

% Plot 1: BER vs SNR
figure;
semilogy(SNRdB_range, BER_results(1,:), 'b-o', ...
         SNRdB_range, BER_results(2,:), 'g-s', ...
         SNRdB_range, BER_results(3,:), 'm-d', ...
         SNRdB_range, BER_adaptive, 'r-*','LineWidth',1.5);
grid on; xlabel('SNR (dB)'); ylabel('BER');
legend('BPSK','QPSK','16QAM','Adaptive');
title('BER vs SNR (Adaptive Modulation)');

% Plot 2: Adaptive Modulation Selection vs SNR
figure;
stairs(SNRdB_range, modulationChoice, 'b','LineWidth',2); grid on;
xlabel('SNR (dB)'); ylabel('Modulation Order (bits/symbol)');
title(['Adaptive Modulation Selection vs SNR (Target BER = ' num2str(targetBER) ')']);
yticks([1 2 4]); yticklabels({'BPSK','QPSK','16QAM'});

% Plot 3: Spectral Efficiency vs SNR
figure;
plot(SNRdB_range, SE_results(1,:),'b-o', ...
     SNRdB_range, SE_results(2,:),'g-s', ...
     SNRdB_range, SE_results(3,:),'m-d', ...
     SNRdB_range, SE_adaptive,'r-*','LineWidth',1.5);
grid on; xlabel('SNR (dB)'); ylabel('Spectral Efficiency (bits/symbol)');
legend('BPSK','QPSK','16QAM','Adaptive');
title('Spectral Efficiency vs SNR');

% Plot 4: Adaptive Modulation BER & SE vs SNR (combined)
figure;
yyaxis left
semilogy(SNRdB_range, BER_adaptive,'b-o','LineWidth',1.5); ylabel('BER');
yyaxis right
plot(SNRdB_range, SE_adaptive,'r-s','LineWidth',1.5); ylabel('Spectral Efficiency');
xlabel('SNR (dB)'); grid on; title('Adaptive Modulation: BER & SE');

% Plot 5: Constellation Diagrams (Fading Effects)
SNRdB_plot = 15;
modulations = {'BPSK','QPSK','16QAM'};

for m = 1:length(modulations)
    symbols = modulate_bits(randi([0 1], 1000,1), modulations{m});

    %Apply different channel effects
    rx_awgn = channel_model(symbols, SNRdB_plot,'AWGN',0);         % AWGN only
    rx_rayleigh = channel_model(symbols, SNRdB_plot,'Rayleigh',0); % Rayleigh fading
    rx_rician = channel_model(symbols, SNRdB_plot,'Rician',ricianK); % Rician fading
    shadowing = 10.^(0.5*randn(size(symbols))/10);
    rx_lognormal = symbols .* shadowing + sqrt(mean(abs(symbols).^2)/10^(SNRdB_plot/10))*(randn(size(symbols))+1j*randn(size(symbols)))/sqrt(2);

    % Subplots
    figure('Name',modulations{m});
    subplot(2,2,1); scatter(real(symbols), imag(symbols),10,'b','filled'); grid on; title('Transmitted Symbols'); xlabel('I'); ylabel('Q'); axis square;
    subplot(2,2,2); scatter(real(rx_awgn), imag(rx_awgn),10,'r','filled'); grid on; title('AWGN'); xlabel('I'); ylabel('Q'); axis square;
    subplot(2,2,3); scatter(real(rx_rayleigh), imag(rx_rayleigh),10,'g','filled'); grid on; title('Rayleigh + AWGN'); xlabel('I'); ylabel('Q'); axis square;
    subplot(2,2,4); scatter(real(rx_rician), imag(rx_rician),10,'m','filled'); grid on; title('Rician + AWGN'); xlabel('I'); ylabel('Q'); axis square;
    sgtitle([modulations{m} ' Constellation: Effect of Fading Channels']);
end

% Plot: BER vs SNR for Fading Channels
fadingTypes = {'AWGN','Rayleigh','Rician'};
BER_fading = zeros(length(schemes), length(SNRdB_range), length(fadingTypes));

fprintf("\nRunning BER vs SNR for each modulation under fading channels...\n");

for f = 1:length(fadingTypes)
    for s = 1:length(schemes)
        scheme = schemes{s};
        for snrIdx = 1:length(SNRdB_range)
            SNRdB = SNRdB_range(snrIdx);
            totalErrors = 0; totalBits = 0;

            for frame = 1:numFrames
                bits = randi([0 1], N, 1);
                symbols = modulate_bits(bits, scheme);
                rxSymbols = channel_model(symbols, SNRdB, fadingTypes{f}, ricianK);
                rxBits = demodulate_symbols(rxSymbols, scheme);
                totalErrors = totalErrors + sum(bits ~= rxBits);
                totalBits = totalBits + length(bits);
            end
            % Store BER
            BER_fading(s, snrIdx, f) = totalErrors / totalBits;
        end
    end
end

% Plot: BER vs SNR for each modulation under fading
for s = 1:length(schemes)
    figure;
    semilogy(SNRdB_range, squeeze(BER_fading(s,:,1)),'b-o','LineWidth',1.5); hold on;
    semilogy(SNRdB_range, squeeze(BER_fading(s,:,2)),'r-s','LineWidth',1.5);
    semilogy(SNRdB_range, squeeze(BER_fading(s,:,3)),'g-d','LineWidth',1.5);
    grid on;
    xlabel('SNR (dB)'); ylabel('BER');
    legend('AWGN','Rayleigh','Rician');
    title(['BER vs SNR for ' schemes{s} ' (With Equalization)']);
end

% Functions
% Channel model: applies AWGN, Rayleigh or Rician fading
function rxSymbols = channel_model(symbols, SNRdB, channelType, ricianK)
    SNR = 10^(SNRdB/10); % Linear SNR
    noise_power = mean(abs(symbols).^2)/SNR;
    switch lower(channelType)
        case 'rayleigh'
            h = (randn(size(symbols)) + 1j*randn(size(symbols)))/sqrt(2); % Rayleigh fading
        case 'rician'
            mean_val = sqrt(ricianK/(ricianK+1));
            std_dev = sqrt(1/(2*(ricianK+1)));
            h = mean_val + std_dev*(randn(size(symbols)) + 1j*randn(size(symbols))); % Rician fading
        case 'awgn'
            h = ones(size(symbols));
        otherwise
            error('Unknown channel type');
    end
    faded_symbols = symbols .* h;
    noise = sqrt(noise_power/2)*(randn(size(symbols))+1j*randn(size(symbols)));
    rxSymbols = faded_symbols + noise;
    if ~strcmpi(channelType,'awgn')
        rxSymbols = rxSymbols ./ h;
    end
end

% Demodulation: converts received symbols back to bits
function bits = demodulate_symbols(symbols, scheme)
    switch scheme
        case 'BPSK'
            bits = real(symbols) > 0;
        case 'QPSK'
            reshaped_bits = [real(symbols)<0, imag(symbols)<0].';
            bits = reshaped_bits(:);
        case '16QAM'
            real_vals = real(symbols)*sqrt(10);
            imag_vals = imag(symbols)*sqrt(10);
            bits = zeros(4,length(symbols));
            bits(1,:) = real_vals < 0;
            bits(2,:) = abs(real_vals) > 2;
            bits(3,:) = imag_vals < 0;
            bits(4,:) = abs(imag_vals) > 2;
            bits = bits(:);
        otherwise
            error('Unknown modulation scheme');
    end
end

% Modulation: converts bits to complex symbols
function symbols = modulate_bits(bits, scheme)
    switch scheme
        case 'BPSK'
            symbols = 2*bits - 1;
        case 'QPSK'
            reshaped_bits = reshape(bits,2,[]);
            symbols = (1/sqrt(2)) * ((1-2*reshaped_bits(1,:)) + 1j*(1-2*reshaped_bits(2,:))).';
        case '16QAM'
            reshaped_bits = reshape(bits,4,[]);
            symbols = zeros(1,size(reshaped_bits,2));
            for i=1:size(reshaped_bits,2)
                real_part = (1-2*reshaped_bits(1,i))*(1+2*reshaped_bits(2,i));
                imag_part = (1-2*reshaped_bits(3,i))*(1+2*reshaped_bits(4,i));
                symbols(i) = (real_part + 1j*imag_part)/sqrt(10);
            end
            symbols = symbols.';
        otherwise
            error('Unknown modulation scheme');
    end
end
