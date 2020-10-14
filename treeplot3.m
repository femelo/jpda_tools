function treeplot3(p,layers,edges,order,c,d)
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
        color = [0.80 0.80 0.80];
        plot(X, Y, '-','Color',color,'LineWidth',2); hold on;
        plot (x, y, 'o','Color',color,'LineWidth',1,'MarkerSize',8,'MarkerFaceColor','w');
    else
        % plot (X, Y, 'g-');
        color = [0.80 0.80 0.80];
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

colors = linspecer(6,'qualitative');
%colors = linspecer(10,'sequential');
%colors = lines(10);
L = cell2mat(layers);
E = zeros(n,n);
for i = 1:n
    for j = 1:n
        if ~isempty(edges{i,j})
            E(i,j) = edges{i,j};
        end
    end
end
layerIds = cell(length(unique(L)),1);
for i = 1:n
    % Annotation of nodes
    tstepid = layers{1,i}+1;
    %cidH = colors(5+tstepid,:);
    cidH = 'k';
    %if layers{1,i} == 0
    %    layerIds{tstepid} = 1:sum(L == tstepid-1);
    %    layerid = layerIds{tstepid}(1);
    %    layerIds{tstepid}(1) = [];
    %    str = sprintf('H_{%d}^{%d}',layerid,tstepid);
    %    fw = 'normal';
    %else
        if isempty(layerIds{tstepid})
            layerIds{tstepid} = 1:sum(L == tstepid-1);
        end
        layerid = layerIds{tstepid}(1);
        layerIds{tstepid}(1) = [];
        if ~transp
            str = sprintf('H_{%d}^{%d}',layerid,tstepid);
        else
            str = sprintf('H_{%d}^{%d}',layerid,tstepid);
        end
        fw = 'normal';
    %end
    text(x(i)+0.005, y(i)+0.01,str, 'Color', cidH, 'FontWeight',fw,'FontSize',15);
    
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
                                str = sprintf('\\theta_{%d}(%d) = %d', layerid, tstepid, edgeid);
                            else
                                str = sprintf('\\theta_{%d}(%d) = %d', layerid, tstepid, edgeid);
                            end
                            fw = 'bold';
                            %fw = 'normal';
                            hp = 'left';
                        case 2
                            % cid = 'r';
                            cid = colors(2,:);
                            if ~transp
                                str = sprintf('\\theta_{%d}(%d) = %d', layerid, tstepid, edgeid);
                            else
                                str = sprintf('\\theta_{%d}(%d) = %d', layerid, tstepid, edgeid);
                            end
                            fw = 'bold';
                            %fw = 'normal';
                            hp = 'left';
                        case 3
                            % cid = 'm';
                            cid = colors(3,:);
                            if ~transp
                                str = sprintf('\\theta_{%d}(%d) = %d', layerid, tstepid, edgeid);
                            else
                                str = sprintf('\\theta_{%d}(%d) = %d', layerid, tstepid, edgeid);
                            end
                            fw = 'bold';
                            %fw = 'normal';
                            hp = 'left';
                        case 4
                            % cid = 'm';
                            cid = colors(4,:);
                            if ~transp
                                str = sprintf('\\theta_{%d}(%d) = %d', layerid, tstepid, edgeid);
                            else
                                str = sprintf('\\theta_{%d}(%d) = %d', layerid, tstepid, edgeid);
                            end
                            fw = 'bold';
                            %fw = 'normal';
                            hp = 'left';
                        case 5
                            % cid = 'm';
                            cid = colors(6,:);
                            if ~transp
                                str = sprintf('\\theta_{%d}(%d) = %d', layerid, tstepid, edgeid);
                            else
                                str = sprintf('\\theta_{%d}(%d) = %d', layerid, tstepid, edgeid);
                            end
                            fw = 'bold';
                            %fw = 'normal';
                            hp = 'left';
                        otherwise
                            %cid = 'k';
                            cid = [0.50 0.50 0.50];
                            %cid = colors(1,:);
                            if ~transp
                                str = sprintf('\\theta_{%d}(%d) = %d', layerid, tstepid, edgeid);
                            else
                                str = sprintf('\\theta_{%d}(%d) = %d', layerid, tstepid, edgeid);
                            end
                            fw = 'normal';
                            hp = 'right';
                    end
                a = 1/4;
                b = 3/4;
                % ft = a + (b-a)*rand(1,1);
                ft = 2/3;
                if x(i) < x(j)
                    xl = x(i)+(x(j)-x(i))*(ft)+0.0030;
                elseif x(i) == x(j)
                    if sum(E(i,:)) > 1
                        xl = x(i)-0.030;
                    else
                        xl = x(i)+0.018;
                    end
                else
                    xl = x(i)-(x(i)-x(j))*(ft)+0.0040;
                end
                yl = y(i)-(y(i)-y(j))*ft;
                text(xl, yl,str,'Color',cid,'FontWeight',fw,'HorizontalAlignment',hp,'FontSize',15);
            end
        end
    end
end



xlabel(['height = ' int2str(h)]);
axis([0 1 0 1]);
set(gca,'Box','off','Visible','off');
