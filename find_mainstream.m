function mainstream = find_mainstream(Q_inlet_win,index_inlet_win,adir_inlet_win)
% traces upstream by following neighbouring cell with highest Q, starting at index_inlet
% returns indices of mainstream cells sorted from lowest to highest index
% can handle inputs w nans but is slower
defdirs

[nrw,ncw]=size(Q_inlet_win);
stream = zeros(nrw,ncw, 'logical');
% [~, outlet]=max(Q_inlet_win{3}(:));
outlet=index_inlet_win;
stream(outlet) = true;
N = nrw*ncw;
pold=0;
Qnbr = zeros(8,1);

% Start at input point and move up
i = outlet;
r = rem(i-1, nrw)+1;
c = fix((i-1)/nrw)+1;
for k = 1:N
    % try to find upstream nbr with highest Q
    Qnbr(:) = 0;
    for d = 1:8
        %Get neighbour in d direction
        rj = r+drow(d);
        cj = c+dcol(d);
        
        % Skip if neighbour is outside domain
        if rj<1,  continue, end;
        if rj>nrw, continue, end;
        if cj<1,  continue, end;
        if cj>ncw, continue, end;
        
        % Skip neighbour if it drains elsewhere
        if adir_inlet_win(rj,cj)~=d, continue, end; % j does not flow to i        
        Qnbr(d) = Q_inlet_win(rj,cj);
    end;
    
    % find max Qnbr
    [Qmax,d] = max(Qnbr); % by default this ignores nans
    if Qmax==0, break; end; % source reached
    
    % set cell upstream as mainstream
    rj = r+drow(d);
    cj = c+dcol(d);
    stream(rj,cj) = true;
    
    % proceed upstream
    r=rj;
    c=cj;
    %fprintf('Progress: %d\n', k);
end;

% figure(2);clf;imagesc(stream);colormap(flipud(gray));axis image;

mainstream=find(stream);

end
