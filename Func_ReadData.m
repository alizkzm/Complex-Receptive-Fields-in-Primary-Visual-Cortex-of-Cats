function output = Func_ReadData(Neuroncode)
    field1 = 'events' ;
    field2 = 'hdr' ;
    value1 = [];
    value2 = [];
    t = 1;
    teststruct = dir(Neuroncode) ;
    output = struct(field1,value1,field2,value2);
    for i = 3 : length(teststruct)
        if(findstr(teststruct(i).name,Neuroncode))
            if(findstr(teststruct(i).name,'msq1D.sa0'))
                 if(findstr(teststruct(i).name,'.sub'))
                     continue ;
                 end
                 if(findstr(teststruct(i).name,'.vecs'))
                     continue ;
                 end
                copyfile ('fget_spk.m' , Neuroncode)
                copyfile ('fget_hdr.m' , Neuroncode )
                cd(teststruct(i).folder);
                output(t).events = (fget_spk(teststruct(i).name))';
                output(t).hdr = fget_hdr(teststruct(i).name);
                t = t+1 ;
            end
            if(findstr(teststruct(i).name,'msq1d.sa0'))
                 if(findstr(teststruct(i).name,'.sub'))
                     continue ;
                 end
                 if(findstr(teststruct(i).name,'.vecs'))
                     continue ;
                 end
                copyfile ('fget_spk.m' , Neuroncode)
                copyfile ('fget_hdr.m' , Neuroncode )
                cd(teststruct(i).folder);
                output(t).events = (fget_spk(teststruct(i).name))';
                output(t).hdr = fget_hdr(teststruct(i).name);
                t = t+1 ;
            end
        end
    end
    %delete 'fget_spk.m'
    %delete 'fget_spk.m'
    cd('C:\Users\Asus\Desktop\Spike_and_Log_Files')
end