
tmp = [linelist.c];
Cent(:,1) = tmp(1:2:end);
Cent(:,2) = tmp(2:2:end);

nodes = [];

for i=1:size(ci_0_to_1,2)
    nodes(i) = find((Cent(:,1)==ci_0_to_1(i).Position(2))&(Cent(:,2)==ci_0_to_1(i).Position(1)));
end
k = i;
for i=1:size(ci_1_to_0,2)
    nodes(i+k) = find((Cent(:,1)==ci_1_to_0(i).Position(2))&(Cent(:,2)==ci_1_to_0(i).Position(1)));
end

