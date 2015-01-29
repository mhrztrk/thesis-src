function save_linelist_info(linelist, filename)

    fid = fopen(filename, 'w+');
    
    for i=1:size(linelist,2)
        fprintf(fid, '%4d - s=[%5.1f,%5.1f] , e=[%5.1f,%5.1f] , clq=[%3d,%3d]\r\n', i, ...
            linelist(i).s(1),linelist(i).s(2),...
            linelist(i).e(1),linelist(i).e(2),...
            linelist(i).clq(1),linelist(i).clq(2));
    end
    
    fclose(fid);

end