function [posds negds] = pick_pos_neg_fromHiC(chrNM)
fn = strcat(chrNM,'.txt');
ds = dataset('File',fn,'HeaderLines',1,'ReadVarNames',false,'ReadObsNames',false);
ds = ds(:,1:7);
ds.Properties.VarNames = {'chr','i','j','x','cor','pos_i','pos_j'};
interact=[ds.i;ds.j];
interact_pos = [ds.pos_i;ds.pos_j];
idd = unique(interact);
freq = [idd, histc(interact, idd)];
[value order] = sort(freq(:,2));
freq_sorted = freq(order, :);
neg_idd = freq_sorted(1:800,1); %negative idx

dssig = ds(ds.cor>=0.4, :);
interact_sig = [dssig.i;dssig.j];
idd_sig = unique(interact_sig);
freq = [idd_sig, histc(interact_sig, idd_sig)];
[value order] = sort(freq(:,2),'descend');
freq_sorted = freq(order, :);
pos_idd = freq_sorted(1:800,1); %positive idx

idd_2_pos = mat2dataset([idd,zeros(length(idd),1)],'VarNames',{'idd','pos'});
for ii = 1:length(idd)
    idd_2_pos.pos(ii) = unique(interact_pos(interact==idd_2_pos.idd(ii)));
end

pos_pos = arrayfun(@(x) idd_2_pos.pos(idd_2_pos.idd == x), pos_idd);
neg_pos = arrayfun(@(x) idd_2_pos.pos(idd_2_pos.idd == x), neg_idd);

posds = mat2dataset([pos_idd, pos_pos],'VarNames',{'idx','pos'});
negds = mat2dataset([neg_idd, neg_pos],'VarNames',{'idx','pos'});
fn = strcat(chrNM,'_pos.txt');
export(posds, 'File', fn, 'Delimiter', ' ');
fn = strcat(chrNM,'_neg.txt');
export(negds, 'File', fn, 'Delimiter', ' ');