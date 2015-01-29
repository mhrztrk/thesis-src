function save_figure(fig, filename)
    tmp = getframe(fig);
    imwrite(tmp.cdata, filename);   
end