function [H,L,S]= HSLPIX(pix)
    MIN=min(pix);
    MAX=max(pix);
    R = pix(1);
    G = pix(2);
    B = pix(3);
    
    L = 0.5*(MAX + MIN);
    if(MAX == MIN)
            H = 0;
            S = 0;
    elseif(R >G && R > B )
            H = mod((60*((G - B)/(MAX - MIN)) + 360),360);
            if(L <= 0.5)
                S = (MAX - MIN)/(2*L);
            else
                S = (MAX - MIN)/(2 - 2*L);
            end
                
    elseif(G > R && G > B)
        H = (60*((B - R)/(MAX - MIN)) + 120);
        if(L <= 0.5)
               S = (MAX - MIN)/(2*L);
         else
                S = (MAX - MIN)/(2 - 2*L);
         end
    else
        H = (60*((R - G)/(MAX - MIN)) + 240);
        if(L <= 0.5)
                S = (MAX - MIN)/(2*L);
         else
                S = (MAX - MIN)/(2 - 2*L);
        end
    end
end
    
    
    