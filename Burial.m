function Burial

class all
clear all
close all
clc

%This code calculates the burial age of a material shielded from cosmic rays.
%It takes into account the elevation of exposure and assumes pre-burial
%exposure was done the method developped by Blard et al.,EPSL, 2019
%The present version of this code is designed for the (10Be-26Al) and
%(10Be-21Ne) pairs.

% Input: Concentrations and uncertainties of 2 radioactive nuclides (10Be
% and 26Al, or 10Be and 21Ne)
% Latitude (°)
% Elevation (m)

addpath(genpath('./functions'));

%Constants
Lambda_Al = log(2)/717000; Lambda_Be = log(2)/1387000; Lambda_Ne = 0;
P_slhl_Al = 27.4; P_slhl_Ne = 17.1; P_slhl_Be = 4.15;

Density = 2.7; Attenuation_length = 160;

Mu = Density / Attenuation_length;

%Initialization
c_names = {'Sample', 'Latitude', 'Altitude', 'X', 's_X', 'Y', 's_Y'};

global Lambda_1 Lambda_2 Lambda_X Lambda_Y P_X P_Y Isotope_X Isotope_Y %input data 

global n Donnees Incertitude_Type Nombre_Colonne_Sortie Sorties Parameter Data_plot %computed

global Sample Lat Elevation X Y D_X D_Y %data


%%%Create the UI 
f = figure('Visible', 'off', 'Position', [200, 200, 900, 500]);
f.Name = 'Burial age calculation';


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

Incertitude_Type = 'Classique'; Nombre_Colonne_Sortie = 5; 
%5 Colonnes : 1)Sample  2)Tburial(Ma) 3)s_Tburial(Ma)  4)Preburial erosion rate 5)1s_Uncertainty

uicontrol('Style', 'pushbutton', 'String', '1 - Calculate burial ages', 'Position', [40, 20, 120, 20], 'Callback', @Calcul);
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
        
        c_names = {'Sample', 'Latitude', 'Present altitude (m)', Isotope_X, ['s_' Isotope_X], Isotope_Y, ['s_' Isotope_Y]};
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
        c_names = {'Sample', 'Latitude', 'Present altitude (m)', Isotope_X, ['s_' Isotope_X], Isotope_Y, ['s_' Isotope_Y]};
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

        X = cell2mat(Donnees.Data(:,4)); 
        D_X = cell2mat(Donnees.Data(:,5)); 
        Y = cell2mat(Donnees.Data(:,6));
        D_Y = cell2mat(Donnees.Data(:,7));

        Parameter = [Lambda_X, P_X, Lambda_Y, P_Y, Mu];
        
        % Check if isotopes are specified
        if(length(Parameter)~=5)
            mydlg = warndlg('Isotopes are not specified !', 'Warning');
            return
        end
       
    end


    function Calcul(~, ~)        
        %Preparing results window
        uicontrol('Style', 'text', 'String', 'Results: ', 'Position' , [400, 450, 200, 20])
        uicontrol('Style', 'pushbutton', 'String', 'Export (.csv)', 'Position', [700, 450, 100, 20], 'Callback', @Export)
        
            Sorties = uitable(f, 'Data', cell(n, Nombre_Colonne_Sortie), 'ColumnName', ...
                {'Sample', 'T_Burial (Ma)', 's_T_Burial (Ma)', 'Preburial_erosion (m/Ma)', '1s'}, ...
                'RowName', [], 'Position', [410, 135, 470, 300], 'ColumnFormat', ...
                {'char', 'numeric', 'numeric', 'numeric', 'numeric'}, ...
                'ColumnEditable', ...
                [false, false, false, false, false]);
            Sorties.Visible = 'on';
        
              
        %Burial age computation and display in a Table
        
           h = waitbar(0,'Computing burial ages...');
        
        for i = drange(1,n)
            
            %close all
            waitbar(i / n)  
            
            if min(Lambda_X, Lambda_Y) == 0
                
                % Case of the 10Be-21Ne pair
                [Tb, s_Tb, Ero, s_Ero] = Burial21Ne_10Be(X(i), D_X(i), Y(i), D_Y(i), Lat(i), Elevation(i), Parameter);
                % export data for plotting later
                Data_plot(i,:) = Tb;
                
            else 
                
                % Case of the 26Al-10Be pair
                [Tb, s_Tb, Ero, s_Ero, T_burial_plot] = Burial26Al_10Be(X(i), D_X(i), Y(i), D_Y(i), Lat(i), Elevation(i), Parameter);
                % export data for plotting later
                Data_plot(i,:) = T_burial_plot;
                                
%                 figure
%                 hist(Tburial_dis,50);
                
            end
            
            Sorties.Data(i,:) = {char(Sample(i)), round(Tb/1e6,2,'significant'), round(s_Tb/1e6,2,'significant'), round(Ero,2,'significant'), round(s_Ero,2,'significant')};
         
        end
        
        close(h);
        
    end


% Banana plot

    function Banane(~,~)
           if min(Lambda_X, Lambda_Y) == 0
            % 10Be-21Ne case
            size(X);
            t_burial = Data_plot(:,1);
            Banana_olives_10_21_burial(Y,D_Y,X,D_X,Lat,Elevation,t_burial,Sample);
        
           else
            % 26Al-10Be case
            size(X);
            t_burial = Data_plot(:,1);
            Banana_olives_26_10_burial(Y,D_Y,X,D_X,Lat,Elevation,t_burial,Sample);
           end
    end

    function Export(~, ~)
            
            S = table(Sorties.Data(:,1), Sorties.Data(:,2), Sorties.Data(:,3), Sorties.Data(:,4), Sorties.Data(:,5));
            S.Properties.VariableNames = {'Sample', 'T_Burial_Ma', 's_T_Burial_Ma','Preburial_erosion_m_per_Ma','One_sigma'};
         
        %Export
        filename = uiputfile({'*.csv'});
        writetable(S, filename, 'Delimiter', '\t')
    end
        
 end

