function Banana_olives_26_10(Al,sAl,Be,sBe,Lat,z,Sample)
% Latitude of normalization // All bananas are drawn for this altitude and
% the measured 10Be concentrations are corrected to be scaled at a 60N latitude :
% 10Be_LatN = 10Be(Sampling_Lat) * Stone(LatN, Sampling_alti) /
% Stone(Sampling_Lat, Sampling_alti)

LatN = 60;

% Production parameters
P26 = 27.4;
P10 = 4.15;
L = 160;
ro = 2.7;
Lambda26 = log(2)/717000;
Lambda10 = log(2)/1387000;


% f is computed using the function: "StoneFactL.m" no VDM assumed for this figure as it doesn'make much sense
gmr = -0.03417;
dtdz = 0.0065;
SLP = 1013.25;

%%%%% Data: 10Be, 1s, 26Al, 1s, latitude, elevation 

% 
% %Sampling elevation
% z = Data(:,6);

% %Sampling latitude 
 Lat = abs(Lat); 


% Normalization of the 10Be and 26Al concentation by the latitude effect, scaling all to Lat 60N

Pressure = SLP .* exp((gmr./dtdz) .* (log(288.15) - log(288.15 - (z.*dtdz))));
Norm = zeros(length(Be),1);

for i = 1:length(Be)
Norm(i) = StoneFactorL(LatN,Pressure(i),SLP)./StoneFactorL(Lat(i),Pressure(i),SLP)';
end

Be = Be.*Norm;
sBe = sBe.*Norm;
Al = Al.*Norm;
sAl = sAl.*Norm;

R = Al./Be;
%sR = ( (sAl./Be).^2 + (sBe.*Al./(Be.^2)).^2).^(1/2);

n = length(R);


%%%%%%%%%%%%%
% 26/10 plot


%Time in ka
t = 0:0.1:100000;

%Erosion in mm/ka or m/Ma
e = 0.000001:0.000001:10;


Z = [0 2000 4000]; %Z = [0:1500:4500];


Pk = SLP .* exp((gmr./dtdz) .* (log(288.15) - log(288.15 - (Z.*dtdz))));

fk = StoneFactorL(LatN,Pk,SLP);

C26_exp = zeros(length(t),length(Z));
C10_exp = zeros(length(t),length(Z));
R_exp = zeros(length(t),length(Z));

% 0-erosion   

for i=1:length(Z)
    
C26_exp(:,i) = fk(i).*(P26./Lambda26)*(1-exp(-Lambda26.*t*1000));
C10_exp(:,i) = fk(i).*(P10./Lambda10)*(1-exp(-Lambda10.*t*1000));
R_exp(:,i) = C26_exp(:,i)./C10_exp(:,i);

end

% steady state erosion

C26_ero = zeros(length(e),length(Z));
C10_ero = zeros(length(e),length(Z));
R_ero = zeros(length(e),length(Z));


for i=1:length(Z)

C26_ero(:,i) = fk(i).*P26./(ro.*e./L+Lambda26);
C10_ero(:,i) = fk(i).*P10./(ro.*e./L+Lambda10);
R_ero(:,i) = C26_ero(:,i)./C10_ero(:,i);

end

  
% Set the Figure
   fignum = figure('Units', 'normal', 'Position', [0.1 0.1 .8 .8]);  hold on %not quite full screen  
   subgroup1 = uipanel('Parent', fignum, 'Units', 'normal', 'Position', [0 0 1 1]);  %top third
   subgroup1_plotbox = uipanel('Parent', subgroup1, 'Units', 'normal', 'Position', [0 0 1 1]);  %plot in top 9/10 of the group
   %subgroup1_controls = uipanel('Parent', subgroup1, 'Units', 'normal', 'Position', [0 0 1 .1]); %control area in bottom 1/10 of the group
   subgroup1_axes = axes('Parent', subgroup1_plotbox, 'Units', 'normal','Position', [0.05 0.1 0.7 0.85]);
   %plot(1:50, rand(1,50), 'Parent', subgroup1_axes);   %throw up some content

   % Create figure
    figure_handle = fignum;
    % create structure of handles
    myhandles = guihandles(figure_handle); 
    % pass all variables to the handles data variable:
     myhandles.Lambda26=Lambda26;
     myhandles.Lambda10=Lambda10;
     myhandles.Be=Be;myhandles.sBe=sBe;
     myhandles.Al=Al;
     myhandles.sAl=sAl;
     myhandles.z = z;
     myhandles.subgroup1_axes = subgroup1_axes;
     % Save the structure
     guidata(figure_handle,myhandles);

   
   % Plot button
   uicontrol('Style', 'pushbutton', 'String', 'Plot selected samples', 'Units', 'normal',...
            'Position', [0.81 0.32 0.12 0.05],'Callback',@plot_list_Callback);
   % listbox   
      subgroup1_listbox = uicontrol('Style','listbox','Parent', subgroup1_plotbox,'string',Sample, ...
       'Units', 'normal','Position', [0.81 0.05 0.15 0.25],'Callback',@listbox_Callback);
      set(subgroup1_listbox,'Max',20,'Min',0);
        listbox_Callback(subgroup1_listbox,[])
   
% Colormap
cc = colormap(othercolor('OrRd9'));
cc = cc(10:64,:);

% color of the curves depending on the altitude
maxz = 4000;
if(max(z)>=4000); maxz = 6000; end
i_z4000 = round(length(cc(:,1))*4000/maxz);
i_z2000 = round(length(cc(:,1))*2000/maxz);

% Drawing of the curves
myhandles.p1 = plot(C10_exp(:,1),R_exp(:,1),'Color',cc(1,:),'LineWidth',3,'Parent', subgroup1_axes); hold on
myhandles.p2 = plot(C10_exp(:,2),R_exp(:,2),'Color',cc(i_z2000,:),'LineWidth',3,'Parent', subgroup1_axes);
myhandles.p3 = plot(C10_exp(:,3),R_exp(:,3),'Color',cc(i_z4000,:),'LineWidth',3,'Parent', subgroup1_axes);

plot(C10_ero(:,1),R_ero(:,1),'Color',cc(1,:),'LineWidth',3,'Parent', subgroup1_axes);
plot(C10_ero(:,2),R_ero(:,2),'Color',cc(i_z2000,:),'LineWidth',3,'Parent', subgroup1_axes);
plot(C10_ero(:,3),R_ero(:,3),'Color',cc(i_z4000,:),'LineWidth',3,'Parent', subgroup1_axes);


%% Plot of data ellipse

 for i = 1: n            
     [ps1,ps2] = data_plot(i,myhandles,subgroup1_axes);
     myhandles.ps1(i)=ps1;
     myhandles.ps2(i)=ps2;     
 end

%% Plot settings
 
% Colorbar
    if(max(z)>=4000)
        caxis([0 6000]);
    else
        caxis([0 4000]);
    end
    h_colorbar = colorbar(subgroup1_axes); 
    ax = gca;
    ax.XScale = 'log';
    ax.YScale = 'log';

 % Axis limits
    min_x = 10^(floor(log10(min(Be))));
    max_x = 10^(ceil(log10(max(C10_exp(:)))));
    max_y = 10^(ceil(log10(max(R_exp(:)))));
    min_y = 10^(floor(log10(min(Al./Be))-1));
    
    xlim(ax,[min_x max_x])
    ylim(ax,[min_y max_y])
    
    % Labels
    xlabel('^{10}Be concentration (at/g) normalized at high latitude >60N');
    ylabel('^{26}Al/^{10}Be');
    ylabel(h_colorbar, 'Present altitude of the site (m)', 'FontSize',11)
    set(gcf, 'Colormap', cc);

    % Color bar
    level = z;
    h_axes = axes('Parent', subgroup1_plotbox, 'position', h_colorbar.Position, 'ylim', h_colorbar.Limits, 'color', 'none', 'visible','off');
    myhandles.p6 = line(h_axes.XLim, level*[1 1], 'color', 'black','parent', h_axes);
    
    % plot legend
    myhandles.leg = legend([myhandles.p1 myhandles.p2 myhandles.p3],{'0 m','2000 m','4000 m'},'Location','southeast');
    title(myhandles.leg,'Paleo-altitude'); 
    
    % Save the handle figure structure
    myhandles.h_axes = h_axes;
    guidata(figure_handle,myhandles);
     
hold off

function [ps1,ps2] = data_plot(i,myhandles,ax)
    % Covariance matrix between Y and X. Caution: only valid for y = b/a vs a
        cor_mat = [[1,0,0,1]; [0,1,0,0]; [0,0,1,0]; [1,0,0,1]];

       Be(i) = myhandles.Be(i);
       sBe(i) = myhandles.sBe(i);
       Al(i) = myhandles.Al(i);
       sAl(i) = myhandles.sAl(i);
       z(i) = myhandles.z(i);
       
    % Interval of confidence
        confidence = 0.68;           
        [u,v] = Ellipse_incertitude_lin([Be(i), sBe(i)], [1, 0],[Al(i), sAl(i)], [Be(i), sBe(i)],cor_mat, confidence ) ;                  
        hold on   
        ps1 = patch(u,v,z(i), 'Parent', ax);
        ps2 = loglog(Be(i),Al(i)/Be(i),'o','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',2, 'Parent', ax);          

function listbox_Callback(hObject, eventdata, handles)
% hObject    handle to input_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
myhandles = guidata(gcbo);
myhandles.index_selected = get(hObject,'Value');
%assignin('base','plot_index_selected',index_selected);
% store data 
guidata(gcbo,myhandles) 

function plot_list_Callback(hObject, eventdata, handles)
% plot selected samples

% get data
myhandles = guidata(gcbo);

% get list of selected samples
list = myhandles.index_selected;

% remove previous curves
delete(myhandles.ps1(:)); delete(myhandles.ps2(:)); 
delete(myhandles.p6(:));
% plot
for i = 1:length(list)
    
    k = list(i);
    
    % altitude on the Colorbar 
    level = myhandles.z(k);
    p6 = line(myhandles.h_axes.XLim, level*[1 1], 'color', 'black','parent', myhandles.h_axes);
    myhandles.p6(i) = p6;
    
    % elipse
     [ps1,ps2] = data_plot(k,myhandles,myhandles.subgroup1_axes);
     myhandles.ps1(i)=ps1;
     myhandles.ps2(i)=ps2;  
end    
     % plot legend
    delete(myhandles.leg);
    myhandles.leg = legend([myhandles.p1 myhandles.p2 myhandles.p3],{'0 m','2000 m','4000 m'},'Location','southeast');
    title(myhandles.leg,'Paleo-altitude');
    
    % save data
    guidata(gcbo,myhandles)

