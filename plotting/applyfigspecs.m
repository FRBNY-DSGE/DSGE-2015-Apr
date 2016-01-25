            
%% Figure Specifications

% Adjust Yaxis
YL = ylabel(Yaxis.axislabel);
set(YL,'FontSize',Yaxis.axislabelsize);
if exist('Yaxis.lim') && ~isempty(Yaxis.lim)
    set(gca,'YLim',Yaxis.lim);
    set(gca,'YTick',(Yaxis.lim(1):Yaxis.freq:Yaxis.lim(2)),'FontSize',Yaxis.ticklabelsize);
end
% 
% if ~isempty(Yaxis.lim)
%     if exist('experiment_flag','var') && experiment_flag==1 && strcmp(plotList{plotnum},'Exp_Forecast')
%         switch V_a 
%             case 1, Yaxis.lim = [-7 16];
%             case 2, Yaxis.lim = [-16 22];
%             case 5, Yaxis.lim = [0 7];
%         end
%     end     
%     set(gca,'YLim',Yaxis.lim);
%     set(gca,'YTick',(Yaxis.lim(1):Yaxis.freq:Yaxis.lim(2)),'FontSize',Yaxis.ticklabelsize);
% end
          
% Adjust Xaxis
set(gca,'XLim',Xaxis.limits);
set(gca,'XTick',Xaxis.freq);
set(gca,'XTickLabel',Xaxis.ticklabels);

% Secondary Yaxis - must be done after Xaxis, but before title to avoid
% resizing
h1 = gca;
h2 = axes('Position',get(gca,'Position'));
set(h2,'YAxisLocation','right','Color','none','XTick',[],'XTickLabel',[]);

% if ~isempty(Yaxis.lim)
%     set(h2,'YLim',Yaxis.lim);
%     set(gca,'YTick',(Yaxis.lim(1):Yaxis.freq:Yaxis.lim(2)),'FontSize',Yaxis.ticklabelsize);
% end

% These lines will align the size, ticks, and fontsize of the two axes.
h(1) = h1;
h(2) = h2;
linkprop(h,{'YLim','YTick','FontSize'});
% for some reason, the YTicks don't synchronize, so this additional line is needed
% if graphs still don't look right, set the rest of the Yaxis properties after the line below
if exist('Yaxis.lim') && ~isempty(Yaxis.lim)
    set(h,'YTick',(Yaxis.lim(1):Yaxis.freq:Yaxis.lim(2)));
end;

% Title
if ~isequal(V_a,Vseq_alt) && ~noTitles
    T = title(Title.name);
    set(T,'FontSize',Title.size);
end

% Misc
if any(strcmp(plotList{plotnum},{'Forecast','Forward Guidance','Forecast Semicond'}))
    set(gca,'Layer','top');
elseif any(strcmp(plotList{plotnum},{'Shock Decomposition'}))
    set(gca,'YGrid','on')
elseif any(strcmp(plotList{plotnum},{'Counterfactual by Variable','Counterfactual by Shock'}))
    grid on;
elseif any(strcmp(plotList,'Shock'))
    grid off
end

