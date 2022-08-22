% Plot a line and points
figure
plot(1:10,1:10,'-o','buttondownfcn',{@Mouse_Callback,'down'});
hold on;
plot(3:12,1:10,'-x');
% Callback function for each point
function Mouse_Callback(hObj,~,action)
persistent curobj xdata ydata ind lag
pos = get(gca,'CurrentPoint');
switch action
    case 'down'
        curobj = hObj;
        xdata = get(hObj,'xdata');
        ydata = get(hObj,'ydata');
        [~,ind] = min(sum((xdata-pos(1)).^2+(ydata-pos(3)).^2,1));
        set(gcf,...
            'WindowButtonMotionFcn',  {@Mouse_Callback,'move'},...
            'WindowButtonUpFcn',      {@Mouse_Callback,'up'});
    case 'move'
        % horizontal move
        if isempty(lag)
            lag = 0;
        end
        lag = lag + pos(2)-xdata(ind);
        disp(lag);
        xdata = xdata + (pos(2) - xdata(ind));
        set(curobj,'xdata',xdata)
    case 'up'
        set(gcf,...
            'WindowButtonMotionFcn',  '',...
            'WindowButtonUpFcn',      '');
end
end