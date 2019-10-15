function [value,isterminal,direction] = dde_stop_events(t,y,YDEL)
% Locate the time when height passes through zero in a decreasing direction
% and stop integration.

isterminal = ones(size(y)); % stop the integration
direction = zeros(size(y));
value = y - ones(size(y)) * 10^10; % detect y-10^8 = 0