function result = runObserverMulti( s, observerType, mapping )
 %# time definition ( s, observerType='RLS', mapping='integral-finite' )
    time = linspace(0.0, (s.Np-1)*s.Tp, s.Np);        
    %# input signal
    %u = asmatrix(10.0*sign(sin(2*pi*(1.0/8.0)*time)));
    u = 10.0*sign(sin(2*pi*(1.0/8.0)*time));
    %# DC motor virtual input
    v = zeros(2,1);
    %# output signal
    y = zeros(1,s.Np);   
    %# plant state    
    x = zeros(2,s.Np);  
    %# load torque
    Tload = 0.01*(ones(1,s.Np)); 
    %# multi-observer    
    mmObserver = MMObserver(s.N, s.M, s.Tp, s.Ao, s.Bo, s.Co, s.Lo, mapping, observerType, 50.0, 1.0, [100.0, 1.0, 0.02]);    
    %# estimated state
    xe = zeros(3,s.Np);   
    %# simulation
    for n = 1:(s.Np-1)  
        %# virtual input - load torque
        v(1,1) = Tload(1,n);
       %# virtual input - voltage
        v(2,1) = u(1,n);
        %# calculate system state
        x(:,n+1) = x(:,n) + s.Tp*(s.A*x(:,n) + s.B*v);
        y(:,n)   = s.C*x(:,n) + 0.05*(rand()-0.5);
        %# call observer
        xe(:,n) = mmObserver.observe(u(:,n), y(:,n));
    end
%     last iteration
    n = s.Np ;        
    y(:,n)   = s.C*x(:,n) + 0.05*(rand()-0.5)
    xe(:,n) = mmObserver.observe(u(:,n), y(:,n))
            
%     result = dict()
%     result('time') = time;
%     result('e1') = xe(0,:)-x(0,:);
%     result('e2') = xe(1,:)-x(1,:);
%     result('e3') = xe(2,:)-Tload;
%     result('x1') = x(0,:);
%     result('x2') = x(1,:);
%     result('x3') = Tload    ;
%     result('xe1') = xe(0,:);
%     result('xe2') = xe(1,:);
%     result('xe3') = xe(2,:);
%     result('u') = u;
    e1 = xe(1,:)-x(1,:);
    e2 = xe(2,:)-x(2,:);
    e3= xe(3,:)-Tload;
    x1= x(1,:);
    x2=x(2,:);
    x3=Tload;
    xe1=xe(1,:);
    xe2=xe(2,:);
    xe3=xe(3,:);
    k ={'time','e1','e2','e3','x1','x2','x3','xe1','xe2','xe3','u'};
    v={time,e1,e2,e3,x1,x2,x3,xe1,xe2,xe3,u};
    result = containers.Map(k,v)
    return
    end