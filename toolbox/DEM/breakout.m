function breakout(cmd)
%BREAKOUT  Pong-like arcade game, a.k.a. Arkanoid/Brickles
%   Breakout begins with five rows of bricks, each row a different color. 
%   The color order from the bottom up is blue, cyan, green, green and red. 
%   Using a single ball, the player must knock down as many bricks as 
%   possible by using the walls and/or the paddle below to ricochet the 
%   ball against the bricks and eliminate them. Blue bricks earn one point 
%   each, cyan bricks earn two points, green bricks earn three points,
%   yellow bricks earn four and the top-level red bricks score five points 
%   each. To add to the challenge, the paddle shrinks to one-half its size 
%   after the ball has broken through the red row and hit the upper wall. 
%   (adapted from Wikipedia)
%
%   Game Controls:
%   Paddle is controlled with mouse.
%   Nine speed levels, changed using num-keys: 1...9.
%   Sound effects are switched on/off using 'S'.
%   Game is paused either with mouse press or 'P'.
%   
%   Example:
%       breakout    % Start Main Breakout Interface
%   Developed by Per-Anders Ekstr�m, 2003-2007 Facilia AB.
%   New pointer control by Joseph Kirk (jdkirk630@gmail.com)
%
%
% https://uk.mathworks.com/matlabcentral/fileexchange/14152-breakout


if ~nargin
    cmd = 'init';
end
if ~(ischar(cmd)||isscalar(cmd))
    return;
end
global paddle bricks ball paddlePos X Y State Sound Speed ...
    hScore hSplash hSound hSpeed
switch cmd
    case 'init' % Initialize the interface
        if ~isempty(findobj('Tag','Breakout'))
            error('Breakout is already running')
        end
        X = 400; Y = 300;
        scrsz = get(0,'ScreenSize');
        % Initialize figure window
        figure('Name','Breakout',...
            'Numbertitle','off',...
            'Menubar','none',...
            'Pointer','custom',...
            'PointerShapeCData',zeros(16)*NaN,...
            'Color','k',...
            'Tag','Breakout',...
            'DoubleBuffer','on',...
            'Colormap',[0 0 0;1 1 1],...
            'WindowButtonDownFcn',sprintf('%s(''NewGame'')',mfilename),...
            'WindowButtonMotionFcn',sprintf('%s(99)',mfilename),...
            'Position',[(scrsz(3)-X)/2 (scrsz(4)-Y)/2 X Y],...
            'CloseRequestFcn',sprintf('%s(''Stop'');closereq;',mfilename),...
            'KeyPressFcn',sprintf('%s(double(get(gcbf,''Currentcharacter'')))',mfilename));
        % Create The axes
        axes('Units','Normalized',...
            'Position', [0 0 1 1],...
            'Visible','off',...
            'DrawMode','fast');
        axis([0,X,0,Y])
        hScore = text(0,Y-Y/20,' Score: 0',...
            'FontName','FixedWidth',...
            'Color','w',...
            'FontWeight','Bold',...
            'FontUnits','normalized');
        hSpeed = text(X/3,Y-Y/20,'Speed \color[rgb]{.4 0 0}||||',...
            'Color',[.3 .3 .3],...
            'FontWeight','Bold',...
            'Visible','on',...
            'HorizontalAlignment','left',...
            'FontUnits','normalized');
        hSound = text(X/7*4,Y-Y/20,'Sound \color[rgb]{.4 0 0}on',...
            'Color',[.3 .3 .3],...
            'FontWeight','Bold',...
            'Visible','on',...
            'HorizontalAlignment','left',...
            'FontUnits','normalized');
        ball = text(X/2,Y/2,'\bullet',...
            'EraseMode','normal',...
            'Color','w',...
            'FontUnits','normalized',...
            'FontSize',.065,...
            'Visible','off',...
            'HorizontalAlignment','Center',...
            'VerticalAlignment','Middle');
        paddlePos = [X/4,0,X/6,Y/20];
        paddle = rectangle('Position',paddlePos,...
            'FaceColor','m',...
            'EraseMode','background');
        bricks = zeros(5,10);
        Sound = true;
        Speed = 35;
        feval(mfilename,'NewBricks')
        hSplash = text(X/2,Y/2,sprintf('BREAKOUT\n(click to start a new game)'),...
            'FontName','FixedWidth',...
            'FontWeight','Bold',...
            'FontUnits','normalized',...
            'HorizontalAlignment','Center',...
            'color','w');
    case 'NewGame'
        delete(findobj(gcbf,'Tag','Lives'))
        set(hSplash,'Visible','off')
        Lives = zeros(3,1);
        for i=1:3
            Lives(i) = text(X-X/17*i,Y-Y/20,'\heartsuit',...
                'Color','r',...
                'FontUnits','normalized',...
                'FontSize',.065,...
                'HorizontalAlignment','Center',...
                'VerticalAlignment','Middle',...
                'Tag','Lives');
        end
        set(gcbf,'WindowButtonDownFcn',sprintf('%s(''PauseGame'')',mfilename))
        set(ball,'Visible','on')
        feval(mfilename,'NewBricks')
        State = 1;
        Score = 0;
        ballX = X/2; ballY = Y/2; ballDX = 5; ballDY = -5;
        Half = 1;
        paddlePos = [X/4,0,X/6,Y/20];
        while(State)
            t0 = clock;
            if State==1
                if ballX<0
                    ballX = 0; 
                    ballDX = -ballDX;
                    if Sound,soundsc(sin(1:100),5000);end
                elseif ballX>X
                    ballX = X; 
                    ballDX = -ballDX;
                    if Sound,soundsc(sin(1:100),5000);end
                end
                YPos = ceil(ballY/Y*20);
                switch YPos
                    case 0 % Ball down the drain
                        if ~isempty(Lives)
                            ballX = X/2;
                            ballY = Y/2;
                            ballDX = 5;
                            ballDY = -5;
                            set(ball,'Position',[ballX ballY])
                            if Sound,soundsc(sin(1:100),1000);end
                            hpos = get(Lives(end),'Position');
                            while(hpos(1)>=0)
                                hpos(1) = hpos(1)-5;
                                set(Lives(end),'Position',hpos)
                                pause(.01)
                            end
                            delete(Lives(end));Lives(end) = [];
                        else
                            set(hSplash,'String',sprintf('GAME OVER\n(click to start a new game)'),'Visible','on')
                            set(ball,'Visible','off')
                            State = 0;
                            if Sound,soundsc(sin(1:100),1000);end
                        end
                    case 1 % Ball in height of paddle
                        if ballX<paddlePos(1)+paddlePos(3) && ballX>paddlePos(1)
                            ballDY = -ballDY;
                            ballDX = (ballX-(paddlePos(1)+paddlePos(3)/2))/3;
                            if Sound,soundsc(sin(1:100),5000);end
                        end
                    case 21 % Ball hit ceiling
                        ballDY = -ballDY;
                        if Sound,soundsc(sin(1:100),5000);end
                        if Half % Make paddle half if first time
                            paddlePos(3) = paddlePos(3)/2;
                            Half = 0;
                        end
                    case {18 17 16 15 14} % Ball in height of bricks
                        YPos = 19-YPos;
                        if ismember(YPos,[1 2 3 4 5])
                            XPos = max(1,ceil(ballX/X*10));
                            if bricks(YPos,XPos)
                                delete(bricks(YPos,XPos))
                                bricks(YPos,XPos) = 0;
                                ballDY = -ballDY;
                                Score = Score+(6-YPos);
                                if Sound,soundsc(sin(1:100),10000);end
                                set(hScore,'String',sprintf(' Score: %i',Score))
                                if isempty(find(bricks,1))
                                    feval(mfilename,'NewBricks')
                                end
                            end
                        end
                end
                % Move ball
                ballX = ballX + ballDX;
                ballY = ballY + ballDY;
                set(ball,'Position',[ballX ballY])
                while etime(clock,t0)<1/Speed
                    drawnow
                end
            else
                pause(.1)
            end
        end
        set(gcbf,'WindowButtonDownFcn',sprintf('%s(''NewGame'')',mfilename))
    case 'PauseGame'
        if State==1
            set(hSplash,'String',sprintf('PAUSE\n(click to start)'),...
                'Visible','on')
            set(ball,'Visible','off')
            State = 2;
        elseif State==2
            set(hSplash,'Visible','off')
            set(ball,'Visible','on')
            State = 1;
        end
    case 'NewBricks'
        delete(findobj(gcbf,'Tag','Bricks'))
        colors = 'rygcb';
        for i=1:5
            for j=0:9
                bricks(i,j+1) = rectangle(...
                    'Position',[j*X/10,Y-Y/20*(2+i),X/10,Y/20],...
                    'Tag','Bricks',...
                    'FaceColor',colors(i));
            end
        end
    case 'Stop'
        State = 0;
    case 99 % Move Paddle
        %- Joseph Kirk (jdkirk630@gmail.com) ------------------------------
        ptr = get(0,'pointerlocation');
        pos = get(gcf,'position');
        set(0,'pointerlocation', [ptr(1) pos(2)+pos(4)/2]);
        %------------------------------------------------------------------
        cp = get(gca,'currentpoint');
        paddlePos(1) = cp(1)-paddlePos(3)/2;
        set(paddle,'Position',paddlePos)
    case 112 % pause
        feval(mfilename,'PauseGame')
    case 115 % sound
        Sound = ~Sound;
        if Sound
        set(hSound,'String','Sound \color[rgb]{.4 0 0}on')
        else
        set(hSound,'String','Sound \color[rgb]{.4 0 0}off')
        end
    case {49 50 51 52 53 54 55 56 57} % Speed
        Speed = 35+(cmd-52)*5;
        set(hSpeed,'String',['Speed \color[rgb]{.4 0 0}' repmat('|',1,cmd-48)])
end
