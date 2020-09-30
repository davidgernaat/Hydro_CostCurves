function basin = fastfindupstream_lim(acc,fdir,drow,dcol, outlet)
%% fast-find catchment upstream of given outlet position
% Davids original file called fastfindupstream_lim,name changed to replace
% orig w faster code
[nr,nc] = size(acc);

[~,order] = sort(acc(:),'descend');

basin = logical(zeros(nr,nc, 'single')); % reset at empty map
basin(outlet) = true; % outlet is lowest grid cell in catchment
N = nr*nc;
firstk = find(order==outlet); % position within 'order' where outlet is.
pold=0;
counter=0;
for k = firstk : N
    counter=counter+1;
    if counter==1000,  continue, end;
    i = order(k);
    %if k==firstk, assert(i==outlet); end;
    if basin(i), continue, end; % happens at outlet location
    r = rem(i-1, nr)+1;
    c = fix((i-1)/nr)+1;
    % find neighbour
    d = fdir(i);
    rj = r+drow(d);
    cj = c+dcol(d);
    if rj<1,  continue, end;
    if rj>nr, continue, end;
    if cj<1,  continue, end;
    if cj>nc, continue, end;
    % copy catchment membership status from downstream nbr
    basin(r,c) = basin(rj,cj);
    % report
    p = fix(100*k/N);
    if p>pold
        %fprintf('Progress: %d%%.\n', p);
        pold=p;
    end;
end
