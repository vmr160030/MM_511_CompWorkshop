clearvars;

which_trials=1:10;
% Load the data.
f = load('ParasolOff_Sp_20190620Ac3.mat');
b_rate = f.analysisParams.binRate;
stimulus_full = f.stimulusParams.stimulus;
response_full = f.stimulusParams.responses;
% Take only the points when the stimulus was on.
valid_pts = f.stimulusParams.preTime + (1:f.stimulusParams.stimTime);
stimulus_full = stimulus_full(:,valid_pts);
response_full = response_full(:,valid_pts);

stimulus_full = stimulus_full(which_trials,:);
response_full = response_full(which_trials,:);

% Let's plot the amplitude spectrum of the stimulus.
figure(100); clf;
subplot(211)
plot(stimulus_full(1,1:500));

subplot(212)
plot(abs(mean(fft(stimulus_full,[],2))));



%% For our purposes, we need to downsample.
bin_rate = 200;
down_samp = b_rate / bin_rate;

bin_edges = linspace(0, size(stimulus_full,2), size(stimulus_full,2)/down_samp+1);

stimulus = zeros(size(stimulus_full,1),size(stimulus_full,2)/down_samp);
response = zeros(size(stimulus_full,1),size(stimulus_full,2)/down_samp);
for ii = 1 : size(stimulus_full,1)
    stimulus(ii,:) = decimate(stimulus_full(ii,:),down_samp);
    response(ii,:) = decimate(response_full(ii,:),down_samp);
end
stimulus = stimulus / down_samp;
response = response / down_samp;
response(response<0)=0;

%% Do correction in the Fourier domain.
f_uncorrected = real(ifft( mean((fft(response,[],2).*conj(fft(stimulus,[],2))),1) ) );

ft = mean((fft(response,[],2).*conj(fft(stimulus,[],2))),1) ...
    ./ mean(fft(stimulus,[],2).*conj(fft(stimulus,[],2)),1) ;

% Take the inverse transform.
filter_fft = real(ifft(ft));

% Define a cutoff frequency.
freqcutoff = 120;
freqcutoff_adjusted = round(freqcutoff/(b_rate/size(stimulus,2))) ; % this adjusts the freq cutoff for the length
ft(:,1+freqcutoff_adjusted:length(stimulus)-freqcutoff_adjusted) = 0 ; 

LinearFilter = real(ifft(ft)) ;

figure(100); clf;
hold on
plot(f_uncorrected(1:100)/norm(f_uncorrected(1:100)));
plot(filter_fft(1:100)/norm(filter_fft(1:100)));
plot(LinearFilter(1:100)/norm(LinearFilter(1:100)));
hold off;


%% Loop through each epoch and compute the covariance and STA matrices.
n_filter_pts = floor(bin_rate * 0.5);
stim_cov = zeros(n_filter_pts,n_filter_pts);
stim_mean = zeros(n_filter_pts,1);
sta = zeros(n_filter_pts,1);
stc = zeros(n_filter_pts,n_filter_pts);
spike_count = 0;
stimulus_count = 0;
for ii = 1 : size(stimulus,1) % loop through each epoch
    % Create the design matrix for the stimulus, taking only valid points
    S = fliplr(toeplitz(stimulus(ii,n_filter_pts:end)', stimulus(ii,n_filter_pts:-1:1)'));
    
    % Stimulus mean.
    stim_mean = stim_mean + mean(S,1)';
    
    % Stimulus covariance.
    stim_cov = stim_cov + S'*S;
    
    % STA
    sta = sta + S' * response(ii,n_filter_pts:end)';
    
    % STC
    stc = S' * (S .* repmat(response(ii,n_filter_pts:end)',1,n_filter_pts));
    
    spike_count = spike_count + sum(response(ii,n_filter_pts:end));
    stimulus_count = stimulus_count + size(S,1);
end

% Divide mean and covariance matrices by stimulus or spike count.
stim_mean = stim_mean / stimulus_count;
stim_cov = stim_cov/(stimulus_count-1) - stim_mean*stim_mean'*stimulus_count/(stimulus_count-1);
sta = sta / spike_count;
stc = stc/(spike_count-1) - sta*sta' * spike_count/(spike_count-1);

%% MLE:

sta_mle = stim_cov \ sta;

%% Ridge regression:
ridge_param = 10; 
sta_ridge = (stim_cov + ridge_param*eye(n_filter_pts)) \ sta;

figure(101); clf;
hold on
plot(sta_mle/norm(sta_mle),'LineWidth',2);
plot(sta_ridge/norm(sta_ridge),'LineWidth',2);
hold off;

