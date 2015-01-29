
figure; imshow(R2, []); hold on;

for i=1:size(CM,1)
    k = convhull(CM{i}(:,1),CM{i}(:,2)); 
    p = plot(CM{i}(k,2),CM{i}(k,1));
    set(p, 'Color', 'red', 'LineWidth', 1.5);
end

scatter(CC(:,2), CC(:,1),'o','b', 'filled');
