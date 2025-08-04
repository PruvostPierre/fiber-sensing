%% Initialization.m - Initialization of the meca test

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Processing of experimental measurements for a coded DAS system
% Authors: P. Pruvost
% Original code by P. Pruvost - 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function TestInProgress = initialization()

    %% 1-Choix du test (CHANGE ME)
    Config=12;
    Test=8;
    FunctionPath = 'D:\01-PHD_TelecomParis\01-TelecomParis_1ère_année_2025-2026\13-Git\fiber-sensing';  % CHANGEME si besoin
    DataBasePath = 'D:\01-PHD_TelecomParis\01-TelecomParis_1ère_année_2025-2026\12-Essais terrain\01_Lannion\05_Data';
    


    %% 2-Variables (CHANGE ME)
    Lambda = 1550;      % Wavelength [nm]
    n = 1.5;            % Refractive index [s.u]
    CphotoElas = 0.78;  % Coefficient photo elastique []
    Lchannel = 0.51;      % Longeur de segment (x_channel) [m]
    
    rowStudied = 278; %(Reflector)
    timeStudied = 1;   %[s]
    timeRecord = 5;    %[s]


    %% 3-Subfunction Path addition
    DataPath = fullfile(DataBasePath, ['Config' num2str(Config)], ['Test' num2str(Test)]);
    fileDAS = fullfile(DataPath, ['DAS_Config' num2str(Config) '_Test' num2str(Test) '.mat']);
    fileIMC = fullfile(DataPath, ['IMC_Config' num2str(Config) '_Test' num2str(Test) '.mat']);

    % Chargement des données DAS
    if exist(fileDAS, 'file')
        DAS = load(fileDAS);
    else
        error(['Fichier DAS non trouvé : ' fileDAS]);
    end

    % Chargement des données IMC
    if exist(fileIMC, 'file')
        IMC = load(fileIMC);
    else
        error(['Fichier IMC non trouvé : ' fileIMC]);
    end
    
    

    %% 4-Structure
    % A- Data
    TestInProgress.DAS = DAS;
    TestInProgress.IMC = IMC;
    TestInProgress.Config = Config;
    TestInProgress.Test = Test;
    TestInProgress.FunctionPath = FunctionPath;
    TestInProgress.DataPath = DataPath;


    % B-Variables

    
    
    TestInProgress.CphotoElas = CphotoElas;
    TestInProgress.Lchannel = Lchannel;
    TestInProgress.Lambda = Lambda;
    TestInProgress.n = n;
    TestInProgress.rowStudied = rowStudied;
    TestInProgress.timeStudied = timeStudied;
    TestInProgress.timeRecord = timeRecord;
   


    %% 5- Modification de structure 
    fields = fieldnames(TestInProgress.DAS.r);
    for i = 1:numel(fields)
        TestInProgress.DAS.(fields{i}) = TestInProgress.DAS.r.(fields{i});
    end

    % Supprimer la structure r
    TestInProgress.DAS = rmfield(TestInProgress.DAS, 'r');

    % Nettoyage
    champsInutiles = {'Acc_Mono_Z_comment', 'Acc_Tri_X_comment', 'Acc_Tri_Y_comment','Acc_Tri_Z_comment','LoadCell50KN_comment'};
    champsExistants = isfield(TestInProgress.IMC, champsInutiles);
    TestInProgress.IMC = rmfield(TestInProgress.IMC, champsInutiles(champsExistants));
end