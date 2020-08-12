%% S Parameter Notebook

% Test code
orig_data = read(rfdata.data,'0806SQ6N0.S2P');
freq = orig_data.Freq;
data = orig_data.S_Parameters(1,1,:);
fit = rationalfit(freq,data);

generateSPICE(fit,'passive.ckt')

