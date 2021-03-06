function ProcessDscML5pFit(dsc,reduced_data,timevect,pLS)
% @fn ProcessDscML5pFit
% @brief Processes measurement descriptor using 5-parameter maximum
%        likelihood fit
% @param dsc The measurement descriptor to process
% @param reduced_data The data vector after adjusting amplitude and time
%                     limits
% @param timevect The sampling times of the values in reduced_data
% @param pLS The result of the LS fit (initial estimators)
% @return none
% @author Tam�s Virosztek, Budapest University of Technology and Economics,
%         Department of Measurement and Infromation Systems,
%         Virosztek.Tamas@mit.bme.hu

screensize = get(0,'ScreenSize');

%Default optimization parameters
MAX_ITER_DEFAULT = 30;
MAX_FUN_EVALS_DEFAULT = 60;
TOL_FUN_DEFAULT = 0;
LAMBDA_ADJUST = 10;
%Constant to set length of pauses
PAUSE_DEFAULT = 0.02; %[sec]

ml_fit_window = figure ('Name','Maximum Likelihood parameter estimation',...
                        'Position',[screensize(3)*0.25 screensize(4)*0.25 screensize(3)*0.5 screensize(4)*0.5],...
                        'NumberTitle','off');

hTextProgressOfEstimation = uicontrol('Style','text',...
                               'Units','normalized',...
                               'Position',[0.1 0.9 0.3 0.05],...
                               'String','Progress of estimation: ',...                               
                               'BackgroundColor',[0.8 0.8 0.8],...
                               'FontWeight','bold');

hTextCosineCoefficient = uicontrol('Style','text',...
                                   'Units','normalized',...
                                   'Position',[0.05 0.8 0.3 0.04],...
                                   'String','Cosine coefficient (A): ',...
                                   'HorizontalAlignment','left',...
                                   'BackgroundColor',[0.8 0.8 0.8]);
                           
hTextCosineCoefficientValue = uicontrol('Style','text',...
                                   'Units','normalized',...
                                   'Position',[0.35 0.8 0.1 0.04],...
                                   'String',sprintf ('%1.5f',0),...
                                   'HorizontalAlignment','right',...
                                   'BackgroundColor',[0.8 0.8 0.8]);

hTextSineCoefficient = uicontrol('Style','text',...
                                   'Units','normalized',...
                                   'Position',[0.05 0.75 0.3 0.04],...
                                   'String','Sine coefficient (B): ',...
                                   'HorizontalAlignment','left',...
                                   'BackgroundColor',[0.8 0.8 0.8]);

hTextSineCoefficientValue = uicontrol('Style','text',...
                                   'Units','normalized',...
                                   'Position',[0.35 0.75 0.1 0.04],...
                                   'String',sprintf ('%1.5f',0),...
                                   'HorizontalAlignment','right',...
                                   'BackgroundColor',[0.8 0.8 0.8]);

hTextDCComponent = uicontrol('Style','text',...
                             'Units','normalized',...
                             'Position',[0.05 0.7 0.3 0.04],...
                             'String','DC component (C): ',...
                             'HorizontalAlignment','left',...
                             'BackgroundColor',[0.8 0.8 0.8]);

hTextDCComponentValue = uicontrol('Style','text',...
                                   'Units','normalized',...
                                   'Position',[0.35 0.7 0.1 0.04],...
                                   'String',sprintf ('%1.5f',0),...
                                   'HorizontalAlignment','right',...
                                   'BackgroundColor',[0.8 0.8 0.8]);

                         
hTextNormalizedFrequency = uicontrol('Style','text',...
                                  'Units','normalized',...
                                  'Position',[0.05 0.65 0.25 0.04],...
                                  'String','Normalized freq. (f/fs): ',...
                                  'HorizontalAlignment','left',...
                                  'BackgroundColor',[0.8 0.8 0.8]);
                              
hTextNormalizedFrequencyValue = uicontrol('Style','text',...
                                   'Units','normalized',...
                                   'Position',[0.3 0.65 0.15 0.04],...
                                   'String',sprintf ('%1.5f',0),...
                                   'HorizontalAlignment','right',...
                                   'BackgroundColor',[0.8 0.8 0.8]);

                         
hTextNoiseStd = uicontrol('Style','text',...
                          'Units','normalized',...
                          'Position',[0.05 0.6 0.25 0.04],...
                          'String','Std. deviation of noise: (sigma)',...
                          'HorizontalAlignment','left',...
                          'BackgroundColor',[0.8 0.8 0.8]);
                      
hTextNoiseStdValue = uicontrol('Style','text',...
                                   'Units','normalized',...
                                   'Position',[0.3 0.6 0.15 0.04],...
                                   'String',sprintf ('%1.5f',0),...
                                   'HorizontalAlignment','right',...
                                   'BackgroundColor',[0.8 0.8 0.8]);

hTextCostFunction = uicontrol('Style','text',...
                          'Units','normalized',...
                          'Position',[0.05 0.55 0.3 0.04],...
                          'String','Value of ML cost function',...
                          'HorizontalAlignment','left',...
                          'BackgroundColor',[0.8 0.8 0.8],...
                          'FontWeight','bold');
                      
hTextCostFunctionValue = uicontrol('Style','text',...
                                   'Units','normalized',...
                                   'Position',[0.35 0.55 0.1 0.04],...
                                   'String','NaN',...
                                   'HorizontalAlignment','right',...
                                   'BackgroundColor',[0.8 0.8 0.8],...
                                   'FontWeight','bold');

hTextFurtherInformation = uicontrol ('Style','text',...
                                     'Units','normalized',...
                                     'Position',[0.1 0.45 0.3 0.05],...
                                     'String','Further Information',...
                                     'BackgroundColor',[0.8 0.8 0.8],...
                                     'FontWeight','bold');
                                 
hTextLambda = uicontrol('Style','text',...
                          'Units','normalized',...
                          'Position',[0.05 0.4 0.2 0.04],...
                          'String','Lambda: ',...
                          'HorizontalAlignment','left',...
                          'BackgroundColor',[0.8 0.8 0.8]);

hTextLambdaValue = uicontrol('Style','text',...
                                   'Units','normalized',...
                                   'Position',[0.25 0.4 0.2 0.04],...
                                   'String',sprintf('%d',0),...
                                   'HorizontalAlignment','right',...
                                   'BackgroundColor',[0.8 0.8 0.8]);

hTextIterationCounter = uicontrol('Style','text',...
                          'Units','normalized',...
                          'Position',[0.05 0.35 0.3 0.04],...
                          'String','Iteration counter: ',...
                          'HorizontalAlignment','left',...
                          'BackgroundColor',[0.8 0.8 0.8]);

hTextIterationCounterValue = uicontrol('Style','text',...
                                   'Units','normalized',...
                                   'Position',[0.35 0.35 0.1 0.04],...
                                   'String',sprintf('%d',0),...
                                   'HorizontalAlignment','right',...
                                   'BackgroundColor',[0.8 0.8 0.8]);                               

hTextFunEvalsCounter = uicontrol('Style','text',...
                          'Units','normalized',...
                          'Position',[0.05 0.3 0.3 0.04],...
                          'String','Function evaluation counter: ',...
                          'HorizontalAlignment','left',...
                          'BackgroundColor',[0.8 0.8 0.8]);

hTextFunEvalsCounterValue = uicontrol('Style','text',...
                                   'Units','normalized',...
                                   'Position',[0.35 0.3 0.1 0.04],...
                                   'String',sprintf('%d',0),...
                                   'HorizontalAlignment','right',...
                                   'BackgroundColor',[0.8 0.8 0.8]);

                               
hTextAdvancedSettings = uicontrol('Style','text',...
                                 'Units','normalized',...
                                 'Position',[0.6 0.9 0.3 0.05],...
                                 'String','Advanced settings',...
                                 'BackgroundColor',[0.8 0.8 0.8],...
                                 'FontWeight','bold');
                             
hTextMaxIter = uicontrol('Style','text',...
                         'Units','normalized',...
                         'Position',[0.55 0.8 0.3 0.04],...
                         'String','Max. number of iterations: ',...
                         'HorizontalAlignment','left',...
                         'BackgroundColor',[0.8 0.8 0.8]);

hEditMaxIter = uicontrol('Style','edit',...
                         'Units','normalized',...
                         'Position', [0.85 0.8 0.1 0.05],...
                         'String','Default',...
                         'Callback',@EditMaxIter_callback);

hTextMaxFunEvals = uicontrol('Style','text',...
                             'Units','normalized',...
                             'Position',[0.55 0.75 0.3 0.04],...
                             'String','Max. number of cost function evaluations: ',...
                             'HorizontalAlignment','left',...                             
                             'BackgroundColor',[0.8 0.8 0.8]);

hEditMaxFunEvals = uicontrol('Style','edit',...
                             'Units','normalized',...
                             'Position', [0.85 0.75 0.1 0.05],...
                             'String','Default',...
                             'Callback',@EditMaxFunEvals_callback);
                         
hTextTolFun = uicontrol('Style','text',...
                        'Units','normalized',...
                        'Position',[0.55 0.7 0.3 0.04],...
                        'String','Termination tolerance on cost function: ',...
                        'HorizontalAlignment','left',...                             
                        'BackgroundColor',[0.8 0.8 0.8]);
                        
hEditTolFun = uicontrol('Style','edit',...
                        'Units','normalized',...
                        'Position', [0.85 0.7 0.1 0.05],...
                        'String','Default',...
                        'Callback',@EditTolFun_callback);
                    

hPushButtonStartPauseResume = uicontrol ('Style','pushbutton',...
                                 'String','Start iteration',...
                                 'Units','normalized',...
                                 'Position',[0.65 0.4 0.2 0.08],...
                                 'Callback',@StartPauseResume_callback);
                         
hPushButtonStop = uicontrol ('Style','pushbutton',...
                                 'String','Stop iteration',...
                                 'Units','normalized',...
                                 'Position',[0.65 0.3 0.2 0.08],...
                                 'Callback',@Stop_callback);
                             
hPushButtonReset = uicontrol ('Style','pushbutton',...
                              'String','Reset',...
                              'Units','normalized',...
                              'Position',[0.65 0.2 0.2 0.08],...
                              'Callback',@Reset_callback);

hPushButtonCompare = uicontrol('Style','pushbutton',...
                               'String','Compare ML vs LS',...
                               'Units','Normalized',...
                               'Position',[0.65 0.1 0.2 0.08],...
                               'Callback',@Compare_callback);                          

hTextInfoLine = uicontrol('Style','text',...
                          'Units','normalized',...
                          'Position',[0.1 0.02 0.8 0.05],...
                          'String','Optimization status: ',...
                          'BackgroundColor','Green');       
                      
%Setting "global" optimization variables:
p = ones(5,1)*NaN;
CF = NaN; PML = NaN; 
grad = ones(5,1)*NaN;
hess = ones(5)*NaN;
l = 0;
INL = zeros(1,2^dsc.NoB-1);
optim_status = '';
termination_reason = '';
if ~isempty(find(reduced_data ~= round(reduced_data),1));
    set(hTextInfoLine,'String','Error: incorrect measurement data');
    set(hTextInfoLine,'BackgroundColor','Red');
end

options.MaxIter = MAX_ITER_DEFAULT;
options.MaxFunEvals = MAX_FUN_EVALS_DEFAULT;
options.TolFun = TOL_FUN_DEFAULT;

iteration_counter = 0;
evaluation_counter = 0;
                     
%Getting initial estimators:
optim_status = 'Initializing';
UpdateDisplay;
M = length(reduced_data);
p(1:2) = pLS(1:2)/(2^dsc.NoB-2);
p(3) = (pLS(3) - 0.5)*1/(2^dsc.NoB - 2); %0.5->0 2^NoB-1.5-> 1
p(4) = pLS(4);
pure_sinewave = p(1)*cos(((1:M).')*p(4)) + p(2)*sin(((1:M).')*p(4)) + p(3);

%Calculating INL:
display_settings.summary_win = 1; %Summary is displayed
display_settings.results_win = 0; %Results are suppressed
display_settings.warning_dialog = 1; %Warning dialogs will be displayed 
INL = ProcessHistogramTest(dsc,display_settings);

%Treating unestimated transition levels: INL linearly reaches 0
lowest_estimated = find(~isnan(INL),1,'first');
highest_estimated = find(~isnan(INL),1,'last');
INL = [linspace(0,INL(lowest_estimated),lowest_estimated).';INL(lowest_estimated+1:highest_estimated-1);linspace(INL(highest_estimated),0,length(INL)-highest_estimated+1).'];
trans_levels = INL2TransLevels(INL);

%Calculating initial estimator of noise
%To get initial noise estimator, 10000 samples of quantized pure sinewave
%are enough
if (length(pure_sinewave) > 1e4)
    quantized_sinewave = QuantizeSignal(pure_sinewave(1:1e4),trans_levels) + 1;
else
    quantized_sinewave = QuantizeSignal(pure_sinewave,trans_levels) + 1;
end
p(5) = rms(quantized_sinewave - reduced_data(1:length(quantized_sinewave)))/(2^dsc.NoB-2);
%Adjusting initial noise variance estimator: in case of arbitrary low noise
%a higher initial estimator is required to be able to perform computations
if (p(5) < eps)
    p(5) = 0.1/(2^dsc.NoB-2); %0.1 LSB
end
optim_status = 'Initialized';
UpdateDisplay;

    function EditMaxIter_callback(source,eventdata)
        if strcmpi(get(source,'String'),'Default')
            options.MaxIter = MAX_ITER_DEFAULT;
        else
            options.MaxIter = str2double(get(source,'String'));
        end
    end

    function EditMaxFunEvals_callback(source,eventdata)
        if strcmpi(get(source,'String'),'Default')
            options.MaxFunEvals = MAX_FUN_EVALS_DEFAULT;
        else
            options.MaxFunEvals = str2double(get(source,'String'));
        end
    end

    function EditTolFun_callback(source,eventdata)
        if strcmpi(get(source,'String'),'Default')
            options.TolFun = TOL_FUN_DEFAULT;
        else
            options.TolFun = str2double(get(source,'String'));
        end
    end

    function StartPauseResume_callback(source,eventdata)
        if strcmpi(optim_status,'Initialized') %Start iteration
            optim_status = 'Running';
            UpdateDisplay;
            [PML,CF,grad,hess] = EvaluateCF(reduced_data,p,dsc.NoB,INL);
            evaluation_counter = evaluation_counter + 1;
            UpdateDisplay;
            l = max(abs(eig(hess))); %initial lambda is the greatest eigenvalue of the Hess matrix
            iteration_counter = 0;
            RunIterations;            
        elseif strcmpi(optim_status,'Running') %Pause iteration
            optim_status = 'Paused';
            UpdateDisplay;
        elseif strcmpi(optim_status,'Paused') % Resume iteration
            optim_status = 'Running';
            UpdateDisplay;
            RunIterations;
        end
    end

    function Stop_callback(source,eventdata)
        optim_status = 'Terminated';
        termination_reason = 'User';
        UpdateDisplay;
    end

    function Reset_callback(source,eventdata)
        %Getting initial estimators:
        p(1:2) = pLS(1:2)/(2^dsc.NoB-2);
        p(3) = (pLS(3) - 0.5)*1/(2^dsc.NoB - 2); %0.5->0 2^NoB-1.5-> 1
        p(4) = pLS(4);
        pure_sinewave = p(1)*cos(((1:M).')*p(4)) + p(2)*sin(((1:M).')*p(4)) + p(3);

        %Calculating initial estimator of noise
        %To get initial noise estimator, 10000 samples of quantized pure sinewave
        %are enough
        if (length(pure_sinewave) > 1e4)
            quantized_sinewave = QuantizeSignal(pure_sinewave(1:1e4),trans_levels) + 1;
        else
            quantized_sinewave = QuantizeSignal(pure_sinewave,trans_levels) + 1;
        end
        p(5) = rms(quantized_sinewave - reduced_data(1:length(quantized_sinewave)))/(2^dsc.NoB-2);
        %Adjusting initial noise variance estimator: in case of arbitrary low noise
        %a higher initial estimator is required to be able to perform computations
        if (p(5) < eps)
            p(5) = 0.1/(2^dsc.NoB-2); %0.1 LSB
        end
        termination_reason = '';
        iteration_counter = 0;
        evaluation_counter = 0;
        CF = NaN;
        optim_status = 'Initialized';
        
        set(hPushButtonStartPauseResume,'String','Start iteration');
        set(hPushButtonStartPauseResume,'Enable','on');
        UpdateDisplay;
    end

    function Compare_callback(source,eventdata)
        CompareLSvsML(dsc,reduced_data,timevect,pLS,p,INL,hess);
    end

    function RunIterations %Core of iterations
             while strcmpi(optim_status,'Running')
                p_next = p - inv(hess + l*eye(5))*grad;
                while (p_next(5) < 0)
                    l  = l*LAMBDA_ADJUST;
                    p_next = p - inv(hess + l *eye(5))*grad;
                end
                [PML_next,CF_next] = EvaluateCF(reduced_data,p_next,dsc.NoB,INL);
                evaluation_counter = evaluation_counter + 1;
                UpdateDisplay;                
                while (CF_next > CF)
                    l = l*LAMBDA_ADJUST;
                    p_next = p - inv(hess + l*eye(5))*grad;
                    while (p_next(5) < 0)
                        l  = l*LAMBDA_ADJUST;
                        p_next = p - inv(hess + l *eye(5))*grad;
                    end
                    [PML_next,CF_next] = EvaluateCF(reduced_data,p_next,dsc.NoB,INL);
                    evaluation_counter = evaluation_counter + 1;
                    UpdateDisplay;                    
                end
                [PML_next,CF_next,grad_next,hess_next] = EvaluateCF(reduced_data,p_next,dsc.NoB,INL);
                evaluation_counter = evaluation_counter + 1;
                UpdateDisplay;                
                l = l*(1/LAMBDA_ADJUST);
                %Checking termination conditions
                if (iteration_counter == options.MaxIter - 1)
                    termination_reason = 'MaxIter';
                    optim_status = 'Terminated';
                elseif (evaluation_counter > options.MaxFunEvals)
                    termination_reason = 'MaxFunEvals';
                    optim_status = 'Terminated';
                elseif (abs(CF - CF_next) < options.TolFun)
                    termination_reason = 'TolFun';
                    optim_status = 'Terminated';
                end
                p = p_next;
                CF = CF_next;
                PML = PML_next;
                grad = grad_next;
                hess = hess_next;
                iteration_counter = iteration_counter + 1;
                UpdateDisplay;                
            end
    end

    function UpdateDisplay
        set(hTextCosineCoefficientValue,'String',sprintf('%1.5f',p(1)));
        set(hTextSineCoefficientValue,'String',sprintf('%1.5f',p(2)));
        set(hTextDCComponentValue,'String',sprintf('%1.5f',p(3)));
        set(hTextNormalizedFrequencyValue,'String',sprintf('%1.3e',p(4)/(2*pi))); %f/fs is displayed
        set(hTextNoiseStdValue,'String',sprintf('%1.3e',p(5)));
        set(hTextCostFunctionValue,'String',sprintf('%1.3e',CF));
        set(hTextLambdaValue,'String',sprintf('%1.3e',l));
        set(hTextIterationCounterValue,'String',sprintf('%d',iteration_counter));
        set(hTextFunEvalsCounterValue,'String',sprintf('%d',evaluation_counter));
        set(hEditMaxIter,'String',sprintf('%d',options.MaxIter));
        set(hEditMaxFunEvals,'String',sprintf('%d',options.MaxFunEvals));
        set(hEditTolFun,'String',sprintf('%1.3e',options.TolFun));
        if (strcmpi(optim_status,'Initializing'))
            set(hTextInfoLine,'String','Calculating initial parameter estimators using 4 parameter LS fit');
            set(hTextInfoLine,'BackgroundColor','Green');
        elseif (strcmpi(optim_status,'Initialized'))
            set(hTextInfoLine,'String','Initial estimators calculated using 4 parameter LS fit');
            set(hTextInfoLine,'BackgroundColor','Green');
            set(hPushButtonStartPauseResume,'String','Start iteration',...
                                            'Enable','on');
            set(hPushButtonStop,'Enable','off');
            set(hPushButtonReset,'Enable','off');
        elseif (strcmpi(optim_status,'Running'))
            set(hTextInfoLine,'String','Performing iterations...');
            set(hTextInfoLine,'BackgroundColor','Green');
            set(hPushButtonStartPauseResume,'String','Pause iteration',...
                                            'Enable','on');
            set(hPushButtonStop,'Enable','on');
            set(hPushButtonReset,'Enable','off');
        elseif (strcmpi(optim_status,'Paused'))
            set(hTextInfoLine,'String','Iteration paused');
            set(hTextInfoLine,'BackgroundColor','Green');
            set(hPushButtonStartPauseResume,'String','Resume iteration',...
                                            'Enable','on');
            set(hPushButtonStop,'Enable','on');
            set(hPushButtonReset,'Enable','on');                                        
        elseif (strcmpi(optim_status,'Terminated'))
            set(hPushButtonStartPauseResume,'String','Start iteration',...
                                            'Enable','off');
            set(hPushButtonStop,'Enable','off');
            set(hPushButtonReset,'Enable','on');
            if (strcmpi(termination_reason,'MaxIter'))
                set(hTextInfoLine,'String','Iteration terminated: reached maximum number of iterations');
                set(hTextInfoLine,'BackgroundColor','Green');
            elseif (strcmpi(termination_reason,'MaxFunEvals'))
                set(hTextInfoLine,'String','Iteration terminated: reached maximum number of cost function evaluations');
                set(hTextInfoLine,'BackgroundColor','Green');
            elseif (strcmpi(termination_reason,'TolFun'))
                set(hTextInfoLine,'String','Iteration terminated: cost function changes less than specified tolerance');
                set(hTextInfoLine,'BackgroundColor','Green');
            elseif (strcmpi(termination_reason,'User'))
                set(hTextInfoLine,'String','Iteration terminated by user');
                set(hTextInfoLine,'BackgroundColor','Green');
            end
        end
        if (iteration_counter ~= 0)
            set (ml_fit_window,'Name',sprintf('Maximum Likelihood Parameter Estimation: step %d of %d',iteration_counter,options.MaxIter));
        end
        pause(PAUSE_DEFAULT);
    end
end