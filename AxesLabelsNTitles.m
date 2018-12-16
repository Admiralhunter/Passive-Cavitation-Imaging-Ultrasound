function    AxesLabelsNTitles(graphicAxes, xLabeled, yLabeled, GraphTitle, FontSize, FontColor, FontWeight)
%
% AxesLabelsNTitles(graphicAxes, xLabeled, yLabeled, GraphTitle, FontSize, FontColor, FontWeight);
%
% Aid for labeling 2D graphics/plots/images
%
% graphicAxes:      is a handle to the graphics to be modified (often gca)
% xLabeled:         is label for x-axis
% yLabeled:         is label for y-axis
% GraphTitle:       is the title of the graph
% FontSize:         is the size of the font to be used
% FontColor:        is the color of the font (e.g. 'black')
% FontWeight:       is the weight of the font, options are: bold, demi, light (e.g. 'black')
% Written by Oliver D. Kripfgans, PhD; Modified by Kevin Haworth, PhD

% FontSize=12;

if nargin == 6
    set(graphicAxes, 'FontWeight','Light','FontSize', FontSize, 'FontName', 'Helvetica','XColor',FontColor,'YColor',FontColor,'Color','none');
else
    if strncmpi(FontWeight,'bold',4)
        set(graphicAxes, 'FontWeight','Bold','FontSize', FontSize, 'FontName', 'Helvetica','XColor',FontColor,'YColor',FontColor,'Color','none');
    elseif strncmpi(FontWeight,'demi',4)
        set(graphicAxes, 'FontWeight','Demi','FontSize', FontSize, 'FontName', 'Helvetica','XColor',FontColor,'YColor',FontColor,'Color','none');
    else
        set(graphicAxes, 'FontWeight','Light','FontSize', FontSize, 'FontName', 'Helvetica','XColor',FontColor,'YColor',FontColor,'Color','none');
    end
end


    set(xlabel(xLabeled),'FontSize', FontSize, 'FontName','Helvetica','Color',FontColor);
    set(ylabel(yLabeled),'FontSize', FontSize, 'FontName','Helvetica','Color',FontColor);
    title(GraphTitle,'FontSize', FontSize, 'FontName','Helvetica','Color',FontColor);
    set(graphicAxes,'LineWidth',1.5)    
    LineMarkhersH = get(graphicAxes,'Children');
    if isprop(LineMarkhersH,'LineWidth');   set(LineMarkhersH,'LineWidth',2.0); end
    if isprop(LineMarkhersH,'MarkerSize');   set(LineMarkhersH,'MarkerSize',10); end

return
