function pts = pts_in_pdf(hdr)

pts = 0;

for e = 1:hdr.header_data.total_epochs
    pts = pts + hdr.epoch_data{e}.pts_in_epoch;
end

