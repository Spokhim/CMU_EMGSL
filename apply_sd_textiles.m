function uni = apply_sd_textiles(uni)

for ii = 1:32:size(uni,1)
    vec = ii:(ii+32-1);
    for ik = 1:size(uni,2)
        [~,tmp_grid] = gradient(reshape(uni(vec,ik),8,4));
        uni(vec,ik) = reshape(tmp_grid,32,1);
    end
end

end