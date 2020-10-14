function treeplot(p,layers,edges,order,c,d)
%TREEPLOT Plot picture of tree.
%   TREEPLOT(p) plots a picture of a tree given a row vector of
%   parent pointers, with p(i) == 0 for a root. 
%
%   TREEPLOT(P,nodeSpec,edgeSpec) allows optional parameters nodeSpec
%   and edgeSpec to set the node or edge color, marker, and linestyle.
%   Use '' to omit one or both.
%
%   Example:
%      treeplot([2 4 2 0 6 4 6])
%   returns a complete binary tree.
%
%   See also ETREE, TREELAYOUT, ETREEPLOT.

%   Copyright 1984-2009 The MathWorks, Inc. 
%   $Revision: 5.12.4.4 $  $Date: 2010/09/02 13:37:05 $

targets = unique(setdiff(cell2mat(layers),0));
if ~exist('order','var')
    order = 1:numel(targets);
end

% Transposed problem
transp = true;

[x,y,h]=treelayout(p);
f = find(p~=0);
pp = p(f);
X = [x(f); x(pp); NaN(size(f))];
Y = [y(f); y(pp); NaN(size(f))];

X = X(:);
Y = Y(:);

if nargin <= 4,
    n = length(p);
    if n < 500,
        % plot (x, y, 'go', X, Y, 'g-');
        color = [0.75 0.75 0.75];
        plot(X, Y, '-','Color',color,'LineWidth',2); hold on;
        plot (x, y, 'o','Color',color,'LineWidth',1,'MarkerSize',8,'MarkerFaceColor','w');
    else
        % plot (X, Y, 'g-');
        color = [0.75 0.75 0.75];
        plot(X, Y, '-','Color',color,'LineWidth',2);
    end;
else
    [~, clen] = size(c);
    if nargin < 3, 
        if clen > 1, 
            d = [c(1:clen-1) '-']; 
        else
            d = 'r-';
        end;
    end;
    [~, dlen] = size(d);
    if clen>0 && dlen>0
        plot (x, y, c, X, Y, d);
    elseif clen>0,
        plot (x, y, c);
    elseif dlen>0,
        plot (X, Y, d);
    else
    end;
end;

colors = linspecer(4,'qualitative');
% colors = linspecer(4,'sequential');
for i = 1:n
    % Annotation of nodes
    layerid = layers{1,i};
    if layerid == 0
        str = '\0';
        fw = 'normal';
    else
        layerid = find(order == layerid)';
        if ~transp
            str = sprintf('M%d',layerid);
        else
            str = sprintf('T%d',layerid);
        end
        fw = 'bold';
    end
    text(x(i)+0.005, y(i)+0.01,str,'FontWeight',fw,'FontSize',16);
    
    % Annotation of edges
    if i <= size(edges,1)
        for j = 2:n
            edgeid = edges{i,j};
            if ~isempty(edgeid)
                    switch edgeid
                        case 1
                            % cid = 'b';
                            cid = colors(1,:);
                            if ~transp
                                str = sprintf('T%d', edgeid);
                            else
                                str = sprintf('%d', edgeid);
                            end
                            fw = 'bold';
                            hp = 'left';
                        case 2
                            % cid = 'r';
                            cid = colors(2,:);
                            if ~transp
                                str = sprintf('T%d', edgeid);
                            else
                                str = sprintf('%d', edgeid);
                            end
                            fw = 'bold';
                            hp = 'left';
                        case 3
                            % cid = 'm';
                            cid = colors(3,:);
                            if ~transp
                                str = sprintf('T%d', edgeid);
                            else
                                str = sprintf('%d', edgeid);
                            end
                            fw = 'bold';
                            hp = 'left';
                        case 4
                            % cid = 'm';
                            cid = colors(4,:);
                            if ~transp
                                str = sprintf('T%d', edgeid);
                            else
                                str = sprintf('%d', edgeid);
                            end
                            fw = 'bold';
                            hp = 'left';
                        otherwise
                            cid = 'k';
                            if ~transp
                                str = 'FA';
                            else
                                str = '0';
                            end
                            fw = 'normal';
                            hp = 'right';
                    end
                a = 1/4;
                b = 3/4;
                % ft = a + (b-a)*rand(1,1);
                ft = 1/2;
                if x(i) < x(j)
                    xl = x(i)+(x(j)-x(i))*ft;
                else
                    xl = x(i)-(x(i)-x(j))*ft;
                end
                yl = y(i)-(y(i)-y(j))*ft;
                text(xl, yl,str,'Color',cid,'FontWeight',fw,'HorizontalAlignment',hp,'FontSize',16);
            end
        end
    end
end



xlabel(['height = ' int2str(h)]);
axis([0 1 0 1]);
set(gca,'Box','off','Visible','off');
