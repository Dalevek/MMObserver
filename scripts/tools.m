function arr = to_plot_array(v, con)
    arr = squeeze(struct2table(v,'AsArray',true));
    if con > 0
        for k = 1:arr.size+1
            if arr(k) > con
                arr(k) = con;
            end
            if arr(k) < - con
                arr(k) = -con;
            end
        end
    end
    return
end

function signalOut = delay(signalIn, memory, delay)
    signalOut = memory(:,1);
    lastIndex = delay - 1;
    for n = 1:lastIndex+1
        memory(:,n) = memory(:,n+1);
    end
    memory(:,lastIndex+1) = singalIn;
    return;
    
end