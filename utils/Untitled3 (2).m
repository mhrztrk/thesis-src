
DrawResultingNetwork(lineSet{i}, linelist{i}, opt_labels{i}, img{i});

%%
tmp = [cursor_img4_missed.Position];
tmp = [tmp(2:2:end);tmp(1:2:end)]';

tmp = [tmp(1:41,:);tmp(44:end,:)];

scatter(tmp(:,2), tmp(:,1), 'r','o','filled','LineWidth', 2);
%%

tmp = [cursor_img4_false.Position];
tmp = [tmp(2:2:end);tmp(1:2:end)]';
scatter(tmp(:,2), tmp(:,1), 'b','o','filled','LineWidth', 2);
%%

tmp = [cursor_img4_missed_pp.Position];
tmp = [tmp(2:2:end);tmp(1:2:end)]';

scatter(tmp(:,2), tmp(:,1), 'r','o','filled','LineWidth', 2);
%%
tmp = [cursor_img1_missed.Position];
tmp = [tmp(2:2:end);tmp(1:2:end)]';

scatter(tmp(:,2), tmp(:,1), 'r','o','filled','LineWidth', 2);


