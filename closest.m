function [i,y] = closest(Data,Value)

% [i,y] = closest(Data,Value)
% Finds the index and value in DATA that is closest to Value
% Data      Data 
% Value     Value

[~,i] = min(abs(Data-Value));
y = Data(i); % Finds first one only!
