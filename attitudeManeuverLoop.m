%% Run attitudeManeuver until hitting stopping condition

%% Setup workspace
clear all

fileName = ['results/' datestr(datetime, 'mmdd-HHMM' )];
diary ([fileName ' Diary.txt'])


%% Iterate
continuation = 1;
parameters.lNorm = 2;

while(parameters.lNorm < 20)
    attitudeManeuver
    save([fileName 'lNorm' num2str(parameters.lNorm) ' Workspace'])
    pause(.1);
end

display('continuation successfully finished')
diary off
