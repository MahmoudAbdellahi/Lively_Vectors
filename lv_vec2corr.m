function stats = lv_vec2corr(X , Y , xlbl, ylbl )
% takes two vectors/matrices X and y and claculates the correlation between them and
% plots the result .. the second dimenstion is the one we loop on

[R,P] = corr(X  ,Y ,'type','spearman');

R = diag(R);
P = diag(P); % diag will get the correponsing-elements' correlations

if length(R)==1 % one pt .. different figure
    % Create figure
    figure1=figure;
    
    % Create axes
    axes1 = axes('Tag','CorrPlot');
    axis off
    % Set the remaining axes properties
    set(axes1,'Color','none','XTick',[0.45 0.5 0.55 0.6 0.65],'YTick',...
        [100 200 300 400]);
    % Create axes
    axes2 = axes('Position',...
        [0.133600917431193 0.206237424547284 0.814220183486242 0.758551307847091],...
        'Tag','PlotMatrixScatterAx');
    hold(axes2,'on');
    
    % Create plot
    plot(X , Y ,'MarkerSize',2,'Marker','o','LineWidth',4,...
        'LineStyle','none',...
        'Color','k','MarkerFaceColor', 'k');
    
    
    % Create ylabel
    ylabel(ylbl, 'Interpreter','none');
    
    % Create xlabel
    xlabel(xlbl, 'Interpreter','none');
    
    box(axes2,'on');
    % Set the remaining axes properties
    set(axes2,'FontSize',16, 'box','off'); % ,'XGrid','on','YGrid','on'
    % Create textbox
    annotation(figure1,'textbox',...
        [0.76053451733778 0.437143566884526 0.296883027522936 0.139230382293763],...
        'Tag','corrCoefs',...
        'String',['R = ' num2str(round(R,2)) , newline,'P = ' num2str(round(P,2))],...
        'FontSize',16,...
        'FitBoxToText','off',...
        'EdgeColor','none');
else 
    xax = 1:length(R);
    plot(xax, R); signi=zeros(size(R));
    mask = P<0.05;
    hold on,
    axes_h = get(gcf,'CurrentAxes');
    signi(mask==1) = axes_h.YLim(2);
    signi(mask==0) = nan;
    % different color: [0.725490212440491 0.898039221763611 0.756862759590149]
    a = area(xax, signi, 'BaseValue',axes_h.YLim(1),'LineStyle','none', 'FaceColor',[0.235294118523598 0.831372559070587 0.0862745121121407]);
    a.FaceAlpha = 0.2;
    title('Spearman correlation Rho');
    axis([axes_h.XLim(1) axes_h.XLim(2) axes_h.YLim(1) axes_h.YLim(2)]);
    xlabel(xlbl, 'Interpreter','none'); ylabel(ylbl, 'Interpreter','none');
end


stats.R = R;
stats.P = P;
end