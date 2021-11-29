function VF =  plotML_RFdist(a, th, d, trgf, ah, VFsz)
% %% load cell info
% utxt = getML_txt('CyTRain_psy.txt');
% a    = utxt.data{strcmp(utxt.name,'a')};
% th   = utxt.data{strcmp(utxt.name,'th')};
% d    = utxt.data{strcmp(utxt.name,'d')};
% trga = utxt.data{strcmp(utxt.name,'trg_a')};
% trgth= utxt.data{strcmp(utxt.name,'trg_dir')};
% 

if nargin<5
    ah = 0;
end

if nargin<6 | isempty(VFsz)
    VFsz = 250; % visual field size
end

VF   = zeros(2*VFsz+1);

for i = 1:length(a)
    r     = floor((d(i)-1)/2);      % force r to be even, not the best way to do it but good for now
    RF    = getRFpatch(r);
    [x,y] = pol2cart(th(i)*pi/180, a(i));
    VF(round(x)-r+VFsz:round(x)+r+VFsz, round(y)-r+VFsz:round(y)+r+VFsz) = ...
        VF(round(x)-r+VFsz:round(x)+r+VFsz, round(y)-r+VFsz:round(y)+r+VFsz)+RF;
    if trgf == 1
        [x,y] = pol2cart((th(i)+180)*pi/180, a(i));
        VF(round(x)-r+VFsz:round(x)+r+VFsz, round(y)-r+VFsz:round(y)+r+VFsz) = ...
            VF(round(x)-r+VFsz:round(x)+r+VFsz, round(y)-r+VFsz:round(y)+r+VFsz)+RF;
    end
end

VF = VF./length(a);

if ah~=0
    if isempty(ah)
        figure
        set(gcf, 'Units', 'inches', 'Position', [0 0 3.8 3.8])
        ah = axes;
        set(ah, 'PlotBoxAspectRatio', [1 1 1], 'Units', 'Normalized',...
                    'Position', [0.6/3.8 0.6/3.8  3/3.8 3/3.8])
    end
    h = imagesc(1-VF');
    if max(max(VF))<0.4
        set(ah, 'CLim', [0.6 1])
    else
        set(ah, 'CLim', [0 1])
    end    
    set(get(h,'Parent'), 'YDir', 'normal',...
              'XTick', 50:100:2*VFsz+1, 'XTickLabel', int2str([-20:10:20]'),...
              'YTick', 50:100:2*VFsz+1, 'YTickLabel', int2str([-20:10:20]'),...
              'FontSize', get(gca, 'FontSize'))
    colormap('hot')
    xl = get(gca, 'XLim');
    yl = get(gca, 'YLim');
    line(mean(xl), mean(yl), 'Marker', '+', 'MarkerSize', get(gca, 'FontSize'), 'Color', 'k', 'LineWidth', 1)
    
    if isempty(ah)
        xlabel('Degree visual angle', 'FontSize', 14)
        ylabel('Degree visual angle', 'FontSize', 14)
    end
%     hold on
%     axes('Units', 'Normalized', 'Position', [0.6/3.8 0.6/3.8  3/3.8 3/3.8], ...
%           'Color', 'none', 'XLim', [-1 1], 'XTick', [], 'YLim', [-1 1], 'YTick', [])
%     line(0, 0, 'Marker', '+', 'MarkerSize', 15, 'Color', 'k', 'LineWidth', 1)
%     hold off
%     x = colorbar;
%     set(x, 'Location', 'west', 'YDir', 'reverse', 'YTickLabel', num2str((1-get(x, 'YTick'))'))
%     text(VFsz+1, VFsz+1, 'x');
%     set(oh, 'FontSize' 25, 'FontWeight', 'bold')
end
            

% subfunction
function p = getRFpatch(rad)
    p = zeros(2*rad+1);
    cx= rad+1; % center x coordinate 
    cy= rad+1; % center y coordinate
    for i = 1:2*rad+1
        for j = 1:2*rad+1
            if (i-cx)^2+(j-cy)^2 <=rad^2 % if pixel = or smaller than the equation of circle, it's in the circle
                p(i,j) = 1;
            end
        end
    end
return
   


