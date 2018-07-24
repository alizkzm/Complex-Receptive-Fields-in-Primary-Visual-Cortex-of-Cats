function output = spike_count_rate(x)
Neurondata = Func_ReadData(x) ;
L = length(Neurondata) ;
Spikecount = 0 ;
for i = 1 : L
    Spikecount = length(Neurondata(i).events) + Spikecount;
end
output = Spikecount / (L*548.7276) ;
end