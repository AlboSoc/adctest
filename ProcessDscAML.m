function ProcessDscAML(dsc,reduced_data,timevect,pLS)
NUM_OF_FOURIER_COEFFS = 15;
MAX_ITER = 1000;
MAX_CF = 5000;

screen_size = get(0,'ScreenSize');
AML_window = figure('Name','Approximate ML (AML) estimation of signal and quantizer parameters in the reduced parameter space',...
    'Position',[screen_size(3)*0.2 screen_size(4)*0.1 screen_size(3)*0.6 screen_size(4)*0.8],...
    'NumberTitle','off');

hTextTitle = uicontrol('Style','text',...
    'Units','normalized',...
    'Position',[0.1 0.92 0.8 0.08],...
    'BackgroundColor',[0.8 0.8 0.8],...
    'String','Approximate ML (AML) estimation of signal and quantizer parameters in the reduced parameter space',...
    'HorizontalAlignment','center',...
    'FontWeight','bold');

% hTextMaxIterations = uicontrol('Style','text',...
%     'Units','normalized',...
%     'Position',[0.1 0.7 0.4 0.08],...
%     'BackgroundColor',[0.8 0.8 0.8],...
%     'String','Maximum number of iterations: ',...
%     'HorizontalAlignment','left');
% 
% hTextMaxIterations_val = uicontrol('Style','text',...
%     'Units','normalized',...
%     'Position',[0.6 0.7 0.2 0.08],...
%     'BackgroundColor',[0.8 0.8 0.8],...
%     'String',sprintf('%d',MAX_ITER),...
%     'HorizontalAlignment','left');
% 
% hTextMaxCFEvals = uicontrol('Style','text',...
%     'Units','normalized',...
%     'Position',[0.1 0.6 0.4 0.08],...
%     'BackgroundColor',[0.8 0.8 0.8],...
%     'String','Maximum number of cost function evaluations: ',...
%     'HorizontalAlignment','left');
% 
% hTextMaxCFEvals_val = uicontrol('Style','text',...
%     'Units','normalized',...
%     'Position',[0.6 0.6 0.2 0.08],...
%     'BackgroundColor',[0.8 0.8 0.8],...
%     'String',sprintf('%d',MAX_CF),...
%     'HorizontalAlignment','left');

hTextResultOfOptimization = uicontrol('Style','text',...
    'Units','normalized',...
    'Position',[0.1 0.87 0.4 0.05],...
    'BackgroundColor',[0.8 0.8 0.8],...
    'String','Result of optimization: ',...
    'HorizontalAlignment','left');
                            
hTextResultOfOptimization_val = uicontrol('Style','text',...
    'Units','normalized',...
    'Position',[0.6 0.87 0.4 0.05],...
    'BackgroundColor',[0.8 0.8 0.8],...
    'String','N/A',...
    'HorizontalAlignment','left');

hTextNumberOfIterations = uicontrol('Style','text',...
    'Units','normalized',...
    'Position',[0.1 0.84 0.4 0.03],...
    'BackgroundColor',[0.8 0.8 0.8],...
    'String','Number of iterations used: ',...
    'HorizontalAlignment','left');
                            
hTextNumberOfIterations_val = uicontrol('Style','text',...
    'Units','normalized',...
    'Position',[0.6 0.84 0.2 0.03],...
    'BackgroundColor',[0.8 0.8 0.8],...
    'String','N/A',...
    'HorizontalAlignment','left');

hTextNumberOfCFEvals = uicontrol('Style','text',...
    'Units','normalized',...
    'Position',[0.1 0.81 0.4 0.03],...
    'BackgroundColor',[0.8 0.8 0.8],...
    'String','Number of CF evaluations used: ',...
    'HorizontalAlignment','left');
                            
hTextNumberOfCFEvals_val = uicontrol('Style','text',...
    'Units','normalized',...
    'Position',[0.6 0.81 0.2 0.03],...
    'BackgroundColor',[0.8 0.8 0.8],...
    'String','N/A',...
    'HorizontalAlignment','left');

hTextInitialCF = uicontrol('Style','text',...
    'Units','normalized',...
    'Position',[0.1 0.78 0.4 0.03],...
    'BackgroundColor',[0.8 0.8 0.8],...
    'String','Initial value of cost function: ',...
    'HorizontalAlignment','left');
                            
hTextInitialCF_val = uicontrol('Style','text',...
    'Units','normalized',...
    'Position',[0.6 0.78 0.2 0.03],...
    'BackgroundColor',[0.8 0.8 0.8],...
    'String','N/A',...
    'HorizontalAlignment','left');

hTextOptimizedCF = uicontrol('Style','text',...
    'Units','normalized',...
    'Position',[0.1 0.75 0.4 0.03],...
    'BackgroundColor',[0.8 0.8 0.8],...
    'String','Value of cost function after optimization: ',...
    'HorizontalAlignment','left');
                            
hTextOptimizedCF_val = uicontrol('Style','text',...
    'Units','normalized',...
    'Position',[0.6 0.75 0.2 0.03],...
    'BackgroundColor',[0.8 0.8 0.8],...
    'String','N/A',...
    'HorizontalAlignment','left');

hTextElapsedTime = uicontrol('Style','text',...
    'Units','normalized',...
    'Position',[0.1 0.72 0.4 0.03],...
    'BackgroundColor',[0.8 0.8 0.8],...
    'String','Elapsed time: ',...
    'HorizontalAlignment','left');
                            
hTextElapsedTime_val = uicontrol('Style','text',...
    'Units','normalized',...
    'Position',[0.6 0.72 0.2 0.03],...
    'BackgroundColor',[0.8 0.8 0.8],...
    'String','N/A s',...
    'HorizontalAlignment','left');

hTextENOB = uicontrol('Style','text',...
    'Units','normalized',...
    'Position',[0.1 0.69 0.4 0.03],...
    'BackgroundColor',[0.8 0.8 0.8],...
    'String','Effective number of bits: ',...
    'HorizontalAlignment','left',...
    'FontWeight','Bold');
                            
hTextENOB_val = uicontrol('Style','text',...
    'Units','normalized',...
    'Position',[0.6 0.69 0.2 0.03],...
    'BackgroundColor',[0.8 0.8 0.8],...
    'String','N/A',...
    'HorizontalAlignment','left',...
    'FontWeight','Bold');

hTextINL = uicontrol('Style','text',...
    'Units','normalized',...
    'Position',[0.1 0.66 0.4 0.03],...
    'BackgroundColor',[0.8 0.8 0.8],...
    'String','Maximal integral nonlinearity: ',...
    'HorizontalAlignment','left',...
    'FontWeight','Bold');
                            
hTextINL_val = uicontrol('Style','text',...
    'Units','normalized',...
    'Position',[0.6 0.66 0.2 0.03],...
    'BackgroundColor',[0.8 0.8 0.8],...
    'String','N/A',...
    'HorizontalAlignment','left',...
    'FontWeight','Bold');

hAxesModTPlot = axes('Position',[0.1 0.4 0.8 0.2]);
hAxesINLPlot = axes('Position',[0.1 0.1 0.8 0.2]);

%Calculating initial estimators for INL values:
display_settings.summary_win = 1; %Summary is displayed
display_settings.results_win = 0; %Results are suppressed
display_settings.warning_dialog = 1; %Warning dialogs will be displayed 
INL = ProcessHistogramTest(dsc,display_settings);

%Treating unestimated transition levels: INL linearly reaches 0
lowest_estimated = find(~isnan(INL),1,'first');
highest_estimated = find(~isnan(INL),1,'last');
INL = [linspace(0,INL(lowest_estimated),lowest_estimated).';INL(lowest_estimated+1:highest_estimated-1);linspace(INL(highest_estimated),0,length(INL)-highest_estimated+1).'];
trans_levels = INL2TransLevels(INL);

M = length(reduced_data);

%Assembling parameter vector:
p0(1:2) = pLS(1:2)/(2^dsc.NoB-2);
p0(3) = (pLS(3) - 0.5)*1/(2^dsc.NoB - 2); %0.5->0 2^NoB-1.5-> 1
p0(4) = pLS(4);
pure_sinewave = p0(1)*cos(((1:M).')*p0(4)) + p0(2)*sin(((1:M).')*p0(4)) + p0(3);

%Calculating initial estimator of noise
%To get initial noise estimator, 10000 samples of quantized pure sinewave
%are enough
if (length(pure_sinewave) > 1e4)
    quantized_sinewave = QuantizeSignal(pure_sinewave(1:1e4),trans_levels) + 1;
else
    quantized_sinewave = QuantizeSignal(pure_sinewave,trans_levels) + 1;
end
p0(5) = rms(quantized_sinewave - reduced_data(1:length(quantized_sinewave)))/(2^dsc.NoB-2);
%Adjusting initial noise variance estimator: in case of arbitrary low noise
%a higher initial estimator is required to be able to perform computations
if (p0(5) < 0.1/(2^dsc.NoB-2)) %0.1 LSB
    p0(5) = 0.1/(2^dsc.NoB-2); %0.1 LSB
end

%Calculating initial Fourier coefficients:
INLFFT = fft(INL);
p0(6) = INLFFT(1);
for k = 1:(NUM_OF_FOURIER_COEFFS-1)/2
    p0(6+2*k-1) = real(INLFFT(k+1));
    p0(6+2*k) = imag(INLFFT(k+1));
end

InitialCF = EvaluateCF_AML(p0,dsc.data,dsc.NoB);

%Setting options of optimization:
options = optimset('Display','off','PlotFcns',@optimplotfval,'MaxIter',MAX_ITER,'MaxFunEvals',MAX_CF);

%Optimization with respect to the signal and quantizer parameters
start_time = tic;
[p,fval,exitflag,output] = fminsearch(@(p)EvaluateCF_AML(p,dsc.data,dsc.NoB), p0, options);
elapsed_time = toc(start_time);

%Calculating ENOB
p(1:2) = p(1:2)*(2^dsc.NoB - 2);
p(3) = p(3)*(2^dsc.NoB - 2) + 0.5;
estimated_sine_wave = p(1)*cos(((1:M).')*p(4)) + p(2)*sin(((1:M).')*p(4)) + p(3);
residuals = reduced_data(timevect) - estimated_sine_wave(timevect);
ENOB = dsc.NoB - log2(sqrt(1/length(residuals)*(residuals.'*residuals))*sqrt(12));
if (ENOB > dsc.NoB)
    ENOB = dsc.NoB;  %saturating ENOB
end

%Calculating INL
INLFFT = zeros(1,2^dsc.NoB-1);
INLFFT(1) = p(6);
for k = 1:(NUM_OF_FOURIER_COEFFS-1)/2
    INLFFT(k+1) = p(6+2*k-1) + 1i*p(6+2*k);
end

for k = 1:NUM_OF_FOURIER_COEFFS
    INLFFT(end-k+1) = conj(INLFFT(k+1));
end

INL = ifft(INLFFT);

if (max(imag(INL)) > 1e-14)
    warning('Calculated INL is incorrect (the imaginary part is larger than quantization noise)');
end

INL = real(INL);    %Discarding imagiary parts owing to quantization noise
INLMax = max(abs(INL));

set(hTextInitialCF_val,'String',sprintf('%e',InitialCF));
set(hTextOptimizedCF_val,'String',sprintf('%e',fval));
set(hTextNumberOfIterations_val,'String',sprintf('%d',output.iterations));
set(hTextNumberOfCFEvals_val,'String',sprintf('%d',output.funcCount));
set(hTextElapsedTime_val,'String',sprintf('%2.3f s',elapsed_time));
set(hTextENOB_val,'String',sprintf('%2.3f bits',ENOB));
set(hTextINL_val,'String',sprintf('%2.3f LSB',INLMax));

if (1 == exitflag)
    set(hTextResultOfOptimization_val,'String','Algorithm converged to minimum of CF');
elseif (0 == exitflag)
    set(hTextResultOfOptimization_val,'String','Algorithm reached maximum number of iterations or CF evaluations');
else
    set(hTextResultOfOptimization_val,'String','Algorithm terminated by user')
end

%Creating MOD T plot
PhaseValues = zeros(1,length(residuals));
for k = 1:length(timevect)
    PhaseValues(k) = p(4)*timevect(k) - 2*pi*floor((p(4)*timevect(k))/(2*pi)); %fractional phase
end
[PhaseValuesSorted,indeces] = sort(PhaseValues);
ResValuesSorted = residuals(indeces);


%Plotting Mod T diagram of residuals
axes(hAxesModTPlot);
plot(PhaseValuesSorted,ResValuesSorted,'b.');
axis([0, 2*pi, -max(abs(residuals)), max(abs(residuals))]);
xlabel('Phase position of sample [Rad]');
ylabel('Fitting residual [LSB]');
title('Mod T plot of residuals');
hold on;
%Plotting reference sine wave:
resolution = 629;
phi_axis = linspace(0,2*pi,resolution);
reference_sinewave = p(1)*cos(((1:resolution).')*2*pi/resolution) + p(2)*sin(((1:resolution).')*2*pi/resolution) + p(3);
reference_sinewave = min(reference_sinewave,2^dsc.NoB*ones(size(reference_sinewave))); %saturating: high threshold
reference_sinewave = max(reference_sinewave,0*ones(size(reference_sinewave)));         %saturating: low threshold
reference_sinewave = reference_sinewave*max(abs(residuals))/(2^(dsc.NoB-1));           %scaling
reference_sinewave = reference_sinewave - mean(reference_sinewave);                    %removing DC
plot(phi_axis,reference_sinewave,'r');
hold off;

%Plotting estimated INL
axes(hAxesINLPlot);
plot(1:2^dsc.NoB-1,INL);
axis([0, 2^dsc.NoB, -INLMax, INLMax]);
xlabel('Code transition levels');
ylabel('INL [LSB]');
title(sprintf('Estimated integral nonlinearity (Fourier approximation, %d real coefficients)',NUM_OF_FOURIER_COEFFS));
end