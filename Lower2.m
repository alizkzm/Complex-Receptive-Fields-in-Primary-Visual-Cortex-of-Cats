%% Part 2 _ Question 2
tests = dir ;
for i = 3:63
    rate = spike_count_rate (tests(i).name) ;
    if(rate<2)
        disp (tests(i).name)
    end
end 