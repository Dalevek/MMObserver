function signalOut = delay(signalIn, memory, delay)
    signalOut = memory(:,1);
    lastIndex = delay - 1;
    for n = 1:lastIndex
        memory(:,n) = memory(:,n+1);
    end
    memory(:,lastIndex) = signalIn;
    return;
    
end