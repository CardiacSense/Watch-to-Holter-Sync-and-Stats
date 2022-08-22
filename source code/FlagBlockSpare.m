function [Flag_spare]=FlagBlockSpare(Flag,SpareSize)
% FlagBlockSpare - inserting spare to flags. 
% 

Flag_spare = Flag;
for k=2:length(Flag)-1
    if(Flag(k)&& ~Flag(k-1))
        if(SpareSize>k)
            Flag_spare(1:k-1)=1;
        else 
            Flag_spare(k-SpareSize:k-1)=1;
        end
    end
    
    if(~Flag(k)&& Flag(k-1))    
        if((k+SpareSize-1)>length(Flag))
            Flag_spare(k:end)=1;
        else
            Flag_spare(k:k+SpareSize-1)=1;
        end
    end
end
    