function out = token_machine(vectors,target_cnt)
%%
np = size(vectors,2);
count = 1;
offset = 0;
while(count<target_cnt+1)
    for i=1:np
        n = size(vectors{i},1);
        for j=1:n
            
            out(count,:) = [vectors{}];
        end
    end
    offset = offset + 1;
    out(count,:) = 
end
end