function Q = Qinput(t)

global Q1 Q2 

if ((0 <= t) && (t <= 0.8))
    Q = Q1 + t * (Q2 - Q1)/ (0.8);
else if ((0.8 <= t) && (t <= 7.2))
        Q = Q2;
    else
        Q = Q1;
    end
end

    
        