load('enzkindata.mat','S','v')
disp(array2table([S' v'],'VariableNames',{'Substrate','V0'}))
enzkin(S,v)