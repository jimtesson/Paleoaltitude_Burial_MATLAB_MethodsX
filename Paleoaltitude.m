function Paleoaltitude

class all
clear all
close all
clc

%This code calculates elevation from pair of radioactive in 
%situ cosmogenic nuclides,  following the method developped by Blard et al.,EPSL, 2019
%The present version of this code is designed for the (10Be-26Al) and
%(10Be-21Ne) pairs.

addpath(genpath('./functions'));

%Constants
Lambda_Al = log(2)/717000; Lambda_Be = log(2)/1387000; Lambda_Ne = 0;
P_slhl_Al = 27.4; P_slhl_Ne = 17.1; P_slhl_Be = 4.15;

Density = 2.7; Attenuation_length = 160;

Mu = Density / Attenuation_length;

%Initialization
c_names = {'Sample', 'Latitude', 'Present altitude', 'Burial age', 'X', 's_X', 'Y', 's_Y'};

global Lambda_1 Lambda_2 Lambda_X Lambda_Y P_X P_Y Isotope_X Isotope_Y %input data 

global n Donnees Incertitude_Type Nombre_Colonne_Sortie Sorties Parameter %computed

global Sample Lat Elevation T X Y D_X D_Y X_cor Y_cor D_X_cor D_Y_cor %data

global StoneFactor StoneFactor_erosion

%%%Create the UI 
f = figure('Visible', 'off', 'Position', [200, 200, 900, 500]);
f.Name = 'Paleoaltitude';


%Nuclides selection
uicontrol('Style', 'text', 'String', 'Nuclides: ', 'Position' , [0, 450, 200, 20], ...
    'HandleVisibility', 'off');
uicontrol('Style', 'popup', 'String', {'', '26Al', '21Ne', '10Be'}, 'FontSize', 10,'Position', [180, 455, 80, 20], 'Callback', @Choix_Isotope_1);
uicontrol('Style', 'popup', 'String', {'', '26Al', '21Ne', '10Be'}, 'FontSize', 10,'Position', [250, 455, 80, 20], 'Callback', @Choix_Isotope_2);

%Sample number
uicontrol('Style', 'pushbutton', 'String', '+', 'Position', [100, 420, 20, 20], 'Callback', @Ajout_ligne);
uicontrol('Style', 'pushbutton', 'String', '-', 'Position', [130, 420, 20, 20], 'Callback', @Suppr_ligne);
uicontrol('Style', 'text', 'String', 'or', 'Position', [160, 420, 30, 20]);
uicontrol('Style', 'pushbutton', 'String', 'Import data (in .csv, .xls or .xlsx)', 'Position', [200, 420, 150, 20], 'Callback', @Import_donnees);

n = 1;
T_in = cell(n,8);
T_in_format = {'char', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric'};
T_in_editable = [true, true, true, true, true, true, true, true];

Donnees = uitable(f, 'Data', T_in, 'ColumnName', c_names, 'RowName', [], 'Position', [10, 95, 380, 300], ...
    'ColumnFormat', T_in_format, 'ColumnEditable', T_in_editable );

Incertitude_Type = 'Classique'; Nombre_Colonne_Sortie = 13;

uicontrol('Style', 'pushbutton', 'String', '1 - Calculate elevations', 'Position', [40, 20, 120, 20], 'Callback', @Calcul);
uicontrol('Style', 'pushbutton', 'String', '2 - Plot', 'Position', [300, 20, 80, 20], 'Callback', @Banane);


%UI display
f.Visible = 'on';

%Function
    function Choix_Isotope()
        Lambda_X = min(Lambda_1, Lambda_2);
        if Lambda_X == Lambda_Be
            Isotope_X = '10Be'; P_X = P_slhl_Be;
        elseif Lambda_X == Lambda_Ne
            Isotope_X = '21Ne'; P_X = P_slhl_Ne;
        elseif Lambda_X == Lambda_Al
            Isotope_X = '26Al'; P_X = P_slhl_Al;
        else
            Isotope_X = 'X'; P_X = 0; 
        end        
        
        Lambda_Y = max(Lambda_1, Lambda_2);
        if Lambda_Y == Lambda_Be
            Isotope_Y = '10Be'; P_Y = P_slhl_Be;
        elseif Lambda_Y == Lambda_Ne
            Isotope_Y = '21Ne'; P_Y = P_slhl_Ne;
        elseif Lambda_Y == Lambda_Al
            Isotope_Y = '26Al'; P_Y = P_slhl_Al;
        else
            Isotope_Y = 'Y'; P_Y = 0;
        end 
        
    end

    function Choix_Isotope_1(source, ~)
        % User choice loading
        str = source.String;
        val = source.Value;
        
        switch str{val};
            case ''
                Lambda_1 = 0;
            case '10Be'
                Lambda_1 = Lambda_Be;
            case '21Ne'
                Lambda_1 = Lambda_Ne;
            case '26Al'
                Lambda_1 = Lambda_Al;         
        end
        
        Choix_Isotope()
        c_names = {'Sample', 'Latitude', 'Present altitude (m)', 'Burial age (Ma)', Isotope_X, ['s_' Isotope_X], Isotope_Y, ['s_' Isotope_Y]};
        Donnees.ColumnName = c_names;
    end

    function Choix_Isotope_2(source, ~)
        % Use choice loading
        str = source.String;
        val = source.Value;
        
        switch str{val};
            case '10Be'
                Lambda_2 = Lambda_Be;
            case '21Ne'
                Lambda_2 = Lambda_Ne;
            case '26Al'
                Lambda_2 = Lambda_Al; 
        end

        Choix_Isotope()
        c_names = {'Sample', 'Latitude', 'Present altitude (m)', 'Burial age (Ma)', Isotope_X, ['s_' Isotope_X], Isotope_Y, ['s_' Isotope_Y]};
        Donnees.ColumnName = c_names;
    end


    function Ajout_ligne(~, ~)
        Donnees.Data(n+1, :) = cell(1,7);
        n = n+1;
    end


    function Suppr_ligne(~, ~)
        if n > 1
            Donnees.Data = Donnees.Data(1:n-1,:);
            n = n-1;
        end
    end


    function Import_donnees(~, ~)

        filename = uigetfile({'*.csv;*.xls;*.xlsx'});
        if filename==0
            return
        end
        if strfind(filename, 'csv')
            Data = readtable(filename);
            Data = table2cell(Data);
            n = length(Data(:,1));
            Donnees.Data = Data(1:n,:);         
        else
            [~, ~, Data] = xlsread(filename);
            n = length(Data(:,1));
            Donnees.Data = Data(2:n,:);
            n = n-1;
        end

        %Data loading in variables
        Sample = Donnees.Data(:,1); 
        Lat = cell2mat(Donnees.Data(:,2)); 
        Elevation = cell2mat(Donnees.Data(:,3)); 
        
        %Conversion of burial age in million years
        T = 1e6.*cell2mat(Donnees.Data(:,4)); 
        X = cell2mat(Donnees.Data(:,5)); 
        D_X = cell2mat(Donnees.Data(:,6)); 
        Y = cell2mat(Donnees.Data(:,7));
        D_Y = cell2mat(Donnees.Data(:,8));

        Parameter = [Lambda_X, P_X, Lambda_Y, P_Y, Mu];
        
        % Check if isotopes are specified
        if(length(Parameter)~=5)
            mydlg = warndlg('Isotopes are not specified !', 'Warning');
            return
        end
        %Correction of concentrations for radioactive decay during burial
        X_cor = X .* exp(Parameter(1).* T);
        Y_cor = Y .* exp(Parameter(3) .* T);
        D_X_cor = D_X.* exp(Parameter(1).* T); 
        D_Y_cor = D_Y.*exp(Parameter(3).*T);
          
    end

   % function Option(~, ~)
   %     option = figure('Position', [400, 400, 400, 400]);
   %     option.Name = 'Options';
        
   % end

    function Calcul(~, ~)        
        %Preparing results window
        uicontrol('Style', 'text', 'String', 'Results: ', 'Position' , [400, 450, 200, 20])
        uicontrol('Style', 'pushbutton', 'String', 'Export (.csv)', 'Position', [700, 450, 100, 20], 'Callback', @Export)
       
            Sorties = uitable(f, 'Data', cell(n, Nombre_Colonne_Sortie), 'ColumnName', ...
                {'Sample', 'H (m)', 'H_min (m)', 'H_max (m)', 'Tint (Ma)', 's_Tint (Ma)', 'H_erosion (m)', 'H_min_Erosion (m)', 'H_max_Erosion (m)', 'Erosion (m/Ma)', 's_Erosion (m/Ma)','Tint (Ma)', 's_Tint (Ma)'}, ...
                'RowName', [], 'Position', [410, 135, 470, 300], 'ColumnFormat', ...
                {'char', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric'}, ...
                'ColumnEditable', ...
                [false, false, false, false, false, false, false, false, false, false, false, false, false]);
            
            
            Sorties.Visible = 'on';
              
        %Elevation computation and display in a Table
        %Stone factor
        
        StoneFactor = zeros(n,1); StoneFactor_erosion = zeros(n,1);
           h = waitbar(0,'Computing altitudes...');
        
        for i = drange(1,n)
            
            %close all
            waitbar(i / n)  
            
            if min(Lambda_X, Lambda_Y) == 0
                
                % Case of the 10Be-21Ne pair 
                % Continuous exposure       
                [H, H_min, H_max, Tint, s_Tint, StoneFactor(i)] = ...
                Alti21Ne_10Be_expo(X_cor(i), D_X_cor(i), Y_cor(i), D_Y_cor(i), Lat(i), Parameter);
                
                % Steady state erosion
                [H_erosion, H_min_erosion, H_max_erosion, Erosion, s_Erosion, Tint_e, s_Tint_e, StoneFactor_erosion(i)] = ...
                Alti21Ne_10Be_erosion(X_cor(i), D_X_cor(i), Y_cor(i), D_Y_cor(i), Lat(i), Parameter);
                
            else 
                % Case of the 26Al-10Be pair
                % Continuous exposure  
                [H, H_min, H_max, Tint, s_Tint, StoneFactor(i)] = ...
                Alti26Al_10Be_expo(X_cor(i), D_X_cor(i), Y_cor(i), D_Y_cor(i), Lat(i), Parameter);
                
                % Steady state erosion
                [H_erosion, H_min_erosion, H_max_erosion, Erosion, s_Erosion, Tint_e, s_Tint_e, StoneFactor_erosion(i)] = ...
                Alti26Al_10Be_erosion(X_cor(i), D_X_cor(i), Y_cor(i), D_Y_cor(i), Lat(i), Parameter);
            end
            
            Sorties.Data(i,:) = {char(Sample(i)), H, H_min, H_max, round(Tint/1e6,3,'significant'), round(s_Tint/1e6,3,'significant'), ...
                H_erosion, H_min_erosion, H_max_erosion, round(Erosion,3,'significant'), round(s_Erosion,2,'significant'), round(Tint_e/1e6,3,'significant'), round(s_Tint_e/1e6,3,'significant')};
        end
        
        close(h);
        
    end


% Banana plot

    function Banane(~,~)
        
           if min(Lambda_X, Lambda_Y) == 0
            % 10Be-21Ne case
            size(X_cor);
            %Banana_olives_10_21(Be,sBe,Ne,sNe,Lat,z)
            Banana_olives_10_21(Y_cor,D_Y_cor,X_cor,D_X_cor,Lat,Elevation,Sample);
        
           else
            % 26Al-10Be case
            %Banana_olives_26_10(Al,sAl,Be,sBe,Lat,z)
            size(X_cor);
            Banana_olives_26_10(Y_cor,D_Y_cor,X_cor,D_X_cor,Lat,Elevation,Sample);
           end
    end

    function Export(~, ~)
          
            S = table(Sorties.Data(:,1), Sorties.Data(:,2), Sorties.Data(:,3), Sorties.Data(:,4), Sorties.Data(:,5), Sorties.Data(:,6), ...
                Sorties.Data(:,7), Sorties.Data(:,8), Sorties.Data(:,9), Sorties.Data(:,10), Sorties.Data(:,11), Sorties.Data(:,12), Sorties.Data(:,13));        
            S.Properties.VariableNames = {'Sample', 'H_m', 'H_min_m', 'H_max_m', 'Tint_Ma', 's_Tint_Ma','H_Erosion_m', 'H_min_Erosion_m', 'H_max_Erosion_m', 'Erosion_m_per_Ma','s_Erosion_m_per_Ma', 'Tint_e_Ma', 's_Tint_e_Ma'};
            
        %Export
        filename = uiputfile({'*.csv'});
        writetable(S, filename, 'Delimiter', '\t')
        
    end
        
 end

