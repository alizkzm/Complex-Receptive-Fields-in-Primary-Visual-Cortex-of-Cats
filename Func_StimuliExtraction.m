function output = Func_StimuliExtraction(events , stimulus_matrix , frequency )
    events = ceil((events * frequency) / 1e04 ) ;
    output=zeros(16,16,length(events));
    for i = 1:length(events)
        if (events (i)>16 && events(i)<32767)
            for j = 0:15
                output(16-j,1:16,i) = stimulus_matrix(events(i)-j,1:16) ;
            end
        end
    end
end