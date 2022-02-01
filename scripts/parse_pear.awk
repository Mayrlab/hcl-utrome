#!/usr/bin/env awk -f

BEGIN { OFS=","; }
{
    if ($0 ~ /^Forward reads file\.+:/)
    {
        match($0, /([^/]+)_R1\.fastq\.gz$/, sample);
        next;
    }

    if ($0 ~ /^Assembled reads \.+:/)
    {
	match($0, /^[^0-9]+([0-9,]+) \/ ([0-9,]+) \(([^%]+)/, assembled);

	## remove commas
	gsub(/,/, "", assembled[1]);
	gsub(/,/, "", assembled[2]);
	next;
    }
    
    if ($0 ~ /^Discarded reads \.+:/)
    {
	match($0, /^[^0-9]+([0-9,]+) \/ ([0-9,]+) \(([^%]+)/, discarded);

	## remove commas
	gsub(/,/, "", discarded[1]);
	gsub(/,/, "", discarded[2]);
	next;
    }

    if ($0 ~ /^Not assembled reads \.+:/)
    {
	match($0, /^[^0-9]+([0-9,]+) \/ ([0-9,]+) \(([^%]+)/, unassembled);

	## remove commas
	gsub(/,/, "", unassembled[1]);
	gsub(/,/, "", unassembled[2]);
	next;
    }

}
END {
    print sample[1], assembled[2], assembled[1], unassembled[1], discarded[1], assembled[3], unassembled[3], discarded[3];
}
