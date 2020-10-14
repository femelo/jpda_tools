function treeplot2(p,layers,edges,c,d)
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

[x,y,h]=treelayout(p);
f = find(p~=0);
pp = p(f);
X = [x(f); x(pp); NaN(size(f))];
Y = [y(f); y(pp); NaN(size(f))];

X = X(:);
Y = Y(:);

if nargin <= 3,
    n = length(p);
    if n < 500,
        plot (x, y, 'go', X, Y, 'g-');
    else
        plot (X, Y, 'g-');
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

for i = 1:n
    % Annotation of nodes
    layerid = layers{1,i};
    if layerid == 0
        str = '\0';
        fw = 'normal';
    else
        str = sprintf('T%d',layerid);
        fw = 'bold';
    end
    text(x(i)+0.01, y(i)+0.01,str,'FontWeight',fw,'FontSize',14);
    
    % Annotation of edges
        for j = 2:n
            edgeid = edges{i,j};
            if ~isempty(edgeid)
                    switch edgeid
                        case 1
                            cid = 'b';
                            str = sprintf('T%d', edgeid);
                            fw = 'bold';
                            hp = 'left';
                        case 2
                            cid = 'r';
                            str = sprintf('T%d', edgeid);
                            fw = 'bold';
                            hp = 'left';
                        case 3
                            cid = 'm';
                            str = sprintf('T%d', edgeid);
                            fw = 'bold';
                            hp = 'left';
                        otherwise
                            cid = 'k';
                            str = 'FA';
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
                text(xl, yl,str,'Color',cid,'FontWeight',fw,'HorizontalAlignment',hp,'FontSize',14);
            end
        end
end



xlabel(['height = ' int2str(h)]);
axis([0 1 0 1]);
set(gca,'Box','off','Visible','off');
