%% Run attitudeManeuver until hitting stopping condition

%% Setup workspace
clear all

fileName = ['results/' datestr(datetime, 'mmdd-HHMM' )];
diary ([fileName ' Diary.txt'])


%% Iterate
continuation = 1;
parameters.lNorm = 2;

while(parameters.lNorm < 20)
    disp(' ')
    attitudeManeuver
    
    save([fileName ' Workspace'])
end

display('continuation successfull finished')
diary off
