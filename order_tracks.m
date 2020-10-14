function [Omega2, order] = order_tracks(Omega, type)
%ORDER_TRACKS Order tracks for efficient hypotheses management
%   Detailed explanation goes here
%
% Coded by:
% Flavio Eler de Melo (flavio.eler@gmail.com)
% University of Liverpool, September, 2013
%

mk = size(Omega,1);
order = (1:mk)';

if mk < 2
    Omega2 = Omega;
    return;
end

Omega2 = Omega;

tm = cell(mk,1);
tset = [];
for j = 1:mk
    cond = (Omega(j,2:end) == 1);
    ind = find(cond);
    tm{j,1} = ind;
    tset = union(tset,tm{j,1});
end
Nt = length(tset);

switch type
    case 'ehm1'
        %% Original ordering algorithm - ordering measurements
        
        tlist = cell(mk,1);
        mlist = cell(mk,1);
        mcomb = cell(4,2);
        
        for j = 1:mk
            mlist{j,1} = j;
        end
        
        comb = nchoosek(1:mk,2);
        mc = zeros(size(comb,1),1);
        mi = zeros(4,1);
        for j = 1:mk-1
            % For each list of tracks, get the union of correspondent gated targets
            for i = 1:mk
                tlist{i,1} = unique(horzcat(tm{mlist{i,1},1}));
            end
            % Compare the intersection in pairs of measurements' lists
            for k = 1:size(comb,1)
                im1 = intersect(tlist{comb(k,1),1},tlist{comb(k,2),1});
                % im1 = union(tlist{comb(k,1),1},tlist{comb(k,2),1});
                if ~isempty(im1)
                    mc(k,1) = length(im1);
                else
                    mc(k,1) = 0;
                end
            end
            % Consider the lists with maximum intersection of associated common gated
            % measurements (maximum intersection)
            [mm1,im1] = max(mc);
            
            if mm1 == 0
            % if mm1 == Nt
                ind = zeros(1,2);
                k1 = 0;
                k = 1;
                while k1 < 2
                    if ~isempty(mlist{k,1});
                        k1 = k1+1;
                        ind(k1) = k;
                    end
                    k = k+1;
                end
                i1 = ind(1,1) == comb(:,1);
                i2 = ind(1,2) == comb(:,2);
                im1 = find(i1 & i2);
            end
            
            % Assemble potential combinations
            mcomb{1,1} = mlist{comb(im1,1),1}; mcomb{1,2} = mlist{comb(im1,2),1};
            mcomb{2,1} = fliplr(mlist{comb(im1,1),1}); mcomb{2,2} = mlist{comb(im1,2),1};
            mcomb{3,1} = mlist{comb(im1,2),1}; mcomb{3,2} = mlist{comb(im1,1),1};
            mcomb{4,1} = fliplr(mlist{comb(im1,2),1}); mcomb{4,2} = mlist{comb(im1,1),1};
            
            % Check for combinations where the new adjoining targets have minimum
            % intersection of gated measurements
            for k = 1:4
               im2 = intersect(tm{mcomb{k,1}(1,end),1},tm{mcomb{k,2}(1,1),1});
               mi(k,1) = length(im2);
            end
            [~,im2] = min(mi);
            
            % Substitute both previous list by the combined one
            mlist{comb(im1,1),1} = unique(horzcat(mcomb{im2,1},mcomb{im2,2}),'stable');
            mlist{comb(im1,2),1} = [];
        end
        ml = [];
        k = 1;
        while isempty(ml)
            ml = mlist{k,1};
            k = k+1;
        end
        order(1:mk,1) = ml;
        
	
    case 'new' % Ordering algorithm proposed by Flavio Eler de Melo - September 2013
        % Disclaimer:
        % This is an enhancement to software that was generated under licence from QinetiQ Limited 
        % (ref QQC4ISR/UOL/SLA/12-2013). The enhancement relates to a patented algorithm, 
        % Patent Reference: 0315349.1. As specified under the agreement, QinetiQ has an irrevocable 
        % non-exclusive worldwide royalty free licence to use this enhancement.

        %% Ordering algorithm 2 - ordering measurements
        
        mlist = zeros(2*mk-1,1);
        
        % Define middle measurement
        % Find the target that can be assigned to the greatest number of
        % measurements
        
        tau = sum(Omega(:,2:end),1);
        tmm = find(tau == max(tau));
        
        % Criterion 1
        % if the target that can be assigned to the greatest number of
        % measurements is not unique, select the measurement with the
        % maximum number of such targets
        sm = zeros(mk,1);
        for j = 1:mk
            sm(j,1) = size(intersect(tm{j,1},tmm),2);
        end
        
        % Criterion 2
        % Check if it is a unique maximum
        [m1,~] = max(sm);
        im1 = m1 == sm;
        while sum(im1) > 1 % not unique
            tau(tmm) = 0;
            tmm = find(tau == max(tau));
            if sum(tau) > 0
                sm = zeros(mk,1);
                for j = 1:mk
                    sm(j,1) = size(intersect(tm{j,1},tmm),2);
                end
                [m,~] = max(sm);
                im = m == sm;
                im2 = im & im1;
                if sum(im2) > 0
                    im1 = im2;
                end
            else
                sm = zeros(mk,1);
                for j = 1:mk
                    sm(j,1) = size(tm{j,1},2);
                end
                [m,~] = max(sm);
                im = m == sm;
                im2 = im & im1;
                if sum(im2) > 0
                    im = find(im2);
                    im2(im(2:end)) = 0;
                    im1 = im2;
                else
                    im = find(im1);
                    im1(im(2:end)) = 0;
                end
            end
        end
        im = find(im1);
        
        % Start with the list with maximum number of assignments
        sm = zeros(mk,1);
        sm2 = zeros(mk,1);
        tset = [];
        for j = 1:mk
            sm(j,1) = size(tm{j,1},2);
            sm2(j,1) = size(tm{j,1},2);
            tset = union(tset,tm{j,1});
        end
        
        % [m1,~] = max(sm);
        
        % Find if it is the only maximum
        % im1 = find(sm == m1);
        % mi = zeros(length(im1),1);
        
        % if length(im1) > 1
        %     for i = 1:length(im1)
        %         im2 = 0;
        %         for j = 1:mk
        %             if j ~= im1(i)
        %                 im = intersect(tm{im1(i),1},tm{j,1});
        %                 im2 = im2 + size(im,2);
        %             end
        %         end
        %         mi(i,1) = im2;
        %     end
        % end
        
        % [~,im2] = max(mi);
        % im = im1(im2);
        
        mlist(mk) = im;
        
        % Find two tracks
        % tstack = zeros(2,1);
        % istack = 0;
        mc = zeros(mk,1);
        exc = im;
        sm2(im) = 0;
        
        % i0 = (mk+mod(mk,2))/2 +1;
        i0 = mk;
        
        % Find the one with greatest intersection of assignment
        for j = 1:mk
            im1 = intersect(tm{im,1},tm{j,1});
            if ~ismember(j,exc)
                if ~isempty(im1)
                    mc(j,1) = length(im1);
                else
                    mc(j,1) = 0;
                end
            end
        end
        [m1,im1] = max(mc);
        % [m1,im1] = min(mc);
        if m1 ~= 0
            mlist(i0-1,1) = im1;
        else
            [~,im1] = max(sm2);
            mlist(i0+1,1) = im1;
        end
        
        mc(:,1) = zeros(size(mc));
        
        j = 2;
        ind = find(mlist ~= 0);
        while j < mk
            iu = ind(1);
            id = ind(end);
            
            % Exceptions list
            exc = mlist(ind,1);
            iud = [iu; id];
            
            if mk-size(ind) >= 2
                for i = 1:2
                    for k = 1:mk
                        im1 = intersect(tm{mlist(iud(i)),1},tm{k,1});
                        if ~ismember(k,exc)
                            if ~isempty(im1)
                                mc(k,1) = length(im1);
                            else
                                mc(k,1) = 0;
                            end
                        end
                    end
                    [m1,im1] = max(mc);
                    % [m1,im1] = min(mc);
                    
                    % Find if it is the only maximum
                    im2 = find(mc == m1);
                    % mi = zeros(length(im1),1);
                    
                    if length(im2) > 1
                        sm2 = sm;
                        sm2(exc) = 0;
                        if i == 1 % Upper
                            [~,im3] = max(sm2(im2));
                        else % Lower
                            % [~,im3] = min(sm(im2));
                            [~,im3] = max(sm2(im2));
                        end
                        im1 = im2(im3);
                        m1 = mc(im1);
                    end
                    
                    % if m1 ~= 0
                        if i == 1
                            mlist(iud(i)-1,1) = im1;
                        else
                            mlist(iud(i)+1,1) = im1;
                        end
                        
                        exc = union(exc, im1);
                    % end
                    mc(:,1) = zeros(size(mc));
                    
                end
            else
                ilast = setdiff(1:mk,mlist(ind));
                for i = 1:2
                    im1 = intersect(tm{mlist(iud(i)),1},tm{ilast,1});
                    if ~isempty(im1)
                        mc(iud(i),1) = length(im1);
                    else
                        mc(iud(i),1) = 0;
                    end
                end
                [m1,im1] = max(mc);
                % [m1,im1] = min(mc);
                if im1 == iud(1) && m1 > 0
                    mlist(iud(1)-1,1) = ilast;
                elseif im1 == iud(2) && m1 > 0
                    mlist(iud(2)+1,1) = ilast;
                else
                    if sm(ilast) < sm(iud(1))
                        mlist(iud(1)-1,1) = ilast;
                    else
                        mlist(iud(2)+1,1) = ilast;
                    end
                end
                mc(:,1) = zeros(size(mc));
            end
            ind = find(mlist ~= 0);
            j = length(ind);
        end
        
        order(1:mk,1) = mlist((mlist ~= 0),1);
        
    otherwise
        order = (1:mk)';
        % error('Unknown ordering type.');
end

% Ordering tracks
Omega2(:,2:end) = Omega(order, 2:end);

end

