
%%%CHANGE THESE PARAMETERS
pxsize = 0.0967821;
DT = 0.0064;
base_dir = '/mnt/data4/SPT_fidelity2024_data/SPT_fidelity2024_data_raw/recordings/TOMM20';
traj_fname = 'tracks_rad=0.375_th=0.8_dist=0.3_gapdist=0.3_framegap=20';
mask_fname = 'C2-210311_COS67_TOM20-Halo_Mito-mNeonG_untreated_3.czi.tif_avg17_crop_Simple Segmentation_bin_openCirc1px_compstrck_1_1015.tif';
out_dir = '/tmp';


%%%CODE STARTS HERE
tiff_info = imfinfo(sprintf('%s/%s', base_dir, mask_fname));
mask = zeros(size(tiff_info, 1), tiff_info(1).Height, tiff_info(1).Width);
for i=1:size(tiff_info, 1)
    mask(i,:,:) = imread(sprintf('%s/%s', base_dir, mask_fname), i);
end

tab_comps = [];
hulls = cell(max(mask(:)));
for i=1:size(mask,1)
    cur = squeeze(mask(i,:,:));
    for j=unique(cur(cur > 0))'
        if isempty(hulls{j})
            hulls{j} = {};
        end

        idxs = find(cur == j);
        [ix, iy] = ind2sub(size(cur), idxs);
        pts = [ix iy] * pxsize;
        tab_comps = [tab_comps; j i mean(pts)];

        hl = struct();
        hl.frame = i;
        hl.poly = pts(convhull(pts), :);
        hulls{j} = [hulls{j}; hl];
    end
end
tr_ids = unique(tab_comps(:,1));


tab = csvread(sprintf('%s/%s', base_dir, traj_fname), 5, 1);
tab = tab(:, [2 8 5 4 8]);
tab(:,2) = tab(:,2) * DT;

tab2 = [];
for i=unique(tab(:,1))'
    tr = tab(tab(:,1) == i, :);
    [no, idxs] = sort(tr(:,5));
    tab2 = [tab2; tr(idxs, :)];
end
tab = tab2;


tab_structs = cell(length(tr_ids), 1);
ii = 1;
for i=tr_ids'
    cur_hulls = hulls{i};
    frames = cellfun(@(x) x.frame, cur_hulls);
    idxs = unique(tab(tab(:,5) >= min(frames) & tab(:,5) <= max(frames), 1));
    for j=idxs'
        tr = tab(tab(:,1) == j, :);
        for k=1:size(tr,1)
            hidx = find(cellfun(@(x) x.frame, cur_hulls) == tr(k, 5) + 1);
            if ~isempty(hidx) && inpolygon(tr(k,3), tr(k,4), cur_hulls{hidx}.poly(:,1), cur_hulls{hidx}.poly(:,2))
                tab_structs{ii} = [tab_structs{ii}; tr];
                break
            end
        end
    end
    ii = ii + 1;
end

for u=1:length(tr_ids)
    selected_comp = tr_ids(u);

    if isempty(tab_structs{u})
        continue
    end

    for selected_tr=unique(tab_structs{u}(:,1))'
        tr = tab_structs{u}(tab_structs{u}(:,1) == selected_tr, :);

        if size(tr, 1) < 5
            continue
        end

        tr_str = tab_comps(tab_comps(:,1) == selected_comp, :);

        if size(tr_str, 1) < 5
            continue
        end

        tr_str(:, 2) = tr_str(:,2) * DT;
        tr_str = tr_str(tr_str(:,2) >= tr(1,2) & tr_str(:,2) < tr(end,2), :);
        tr = tr(tr(:,2) >= tr_str(1,2) & tr(:,2) <= tr_str(end,2), :);
        tr_str = tr_str(tr_str(:,2) >= tr(1,2) & tr_str(:,2) < tr(end,2), :);
        
        idxs = [];
        j = 1;
        for i=1:size(tr,1)
            while j <= size(tr_str,1) && tr(i,2) ~= tr_str(j,2)
                j = j + 1;
            end
            if j < size(tr_str,1) && tr(i,2) == tr_str(j,2)
                idxs = [idxs; i j];
            end
        end

        figure
        hold on
        plot(tr_str(:,3), tr_str(:,4), 'k')
        plot(tr(:,3), tr(:,4), 'b')
        hold off
        daspect([1 1 1])
        legend({'Structure', 'Particle'})
        pause(1)
        print(sprintf('%s/trajStruct_vs_traj_%d_%d.svg', out_dir, selected_comp, selected_tr), '-dsvg')

        tr2 = [];
        for i=1:(size(tr,1)-2)
            tmp = tr(i,:);
            tmp(3:4) = tmp(3:4) - tr_str(idxs(i, 2), 3:4);
            tr2 = [tr2; tmp];
        end

        figure
        hold on
        plot(tr(:,3) - tr_str(idxs(1,2), 3) , tr(:,4) - tr_str(idxs(1,2), 4), 'k')
        plot(tr2(:,3), tr2(:,4), 'r')
        hold off
        legend({'Raw', 'Decoupled'})
        daspect([1 1 1])
        pause(1)
        print(sprintf('%s/trajs_corr_%d_%d.svg', out_dir, selected_comp, selected_tr), '-dsvg')
    end
end