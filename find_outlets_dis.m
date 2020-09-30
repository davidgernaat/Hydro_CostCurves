function [ro,co] = find_outlets_dis(flowdist,acc,do)

% Find outlets based on hydrological distance from main outlet
acc(isnan(acc))=1;
ma = mean(acc(:));
channels = acc>ma; % only cells above mean flow accumulation area


stations = find(channels & rem(flowdist, do)==0); % select pixels at equal distance intervals AND within channel network
[ro,co] = ind2sub(size(acc), stations);

% figure(1);imagesc(X); hold on; 
% plot(co,ro,'w.', 'markersize',20); 
% hold off;
end