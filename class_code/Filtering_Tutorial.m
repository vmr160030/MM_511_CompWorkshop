
% Define a Gaussian filter
gauss = @(p,x)((p(2)*sqrt(2*pi))^-1 * exp(-0.5*((x+p(1)).^2/(2*p(2)^2))));

% Generate a random signal.
noise_stream = RandStream('mt19937ar','Seed', 1);
s = noise_stream.randn(1,10000);

% Define a time domain for the filter.
t=0:100;
% Create the filter.
filter_mean=0.0;
filter_std = 10;
f = gauss([filter_mean,filter_std],t);


% plot(t,gauss([0,20],t))


%% Task 1: Perform filtering with nested for loops.
% Reflect the filter in the time domain to do this way...
f = fliplr(f);

f1 = zeros(size(s));
for ii = length(f)+1 : length(s)
    for jj = 1 : length(f)
        f1(ii) = f1(ii) + f(jj)*s(ii-length(f)+jj);
    end
end

figure(100); clf;
hold on
plot(s);
plot(f1);
hold off;

%% Task 2: Perform filtering with matrix multiplication.
% Define a time domain for the filter.
t=0:100;
% Create the filter.
f = gauss([filter_mean,filter_std],t);
% Reflect the filter in the time domain to do this way...
f = fliplr(f);

S = zeros(length(s), length(f));
for ii = length(f) : length(s)
    S(ii,:) = s(ii-length(f)+1:ii);
end

% A more elegant way using a toeplitz matrix.
S2 = zeros(size(S));
S2(length(f):end,:) = fliplr(toeplitz(s(length(f):end), s(length(f):-1:1)));

% Take the matrix product to perform convolution.
f2 = S * f(:);

%% Task 3: Perform filtering in the Fourier domain.

% Define a time domain for the filter.
t=0:length(s)-1;
% Create the filter.
f = gauss([filter_mean,filter_std],t);
f3 = real(ifft( fft(f) .* fft(s) ));

% OK. Let's do this right.
pad_pts=100;
f4 = real(ifft( fft([f,zeros(1,pad_pts)]) .* fft([s,zeros(1,pad_pts)]) ));
f4 = f4(1:length(s));


save('conv_exercise.mat','s','f','f1','f2','f3','f4');


