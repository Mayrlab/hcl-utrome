#!/usr/bin/env awk -f

BEGIN { OFS=","; }
{
    if ($0 ~ /^# whitelist[^:]+:/)
    {
        match($0, /([^/]+)\.whitelist\.txt$/, sample);
        next;
    }

    if ($0 ~ /INFO Input Reads:/)
    {
	match($0, /([0-9]+)$/, total_in);
	next;
    }

    if ($0 ~ /INFO Reads output:/)
    {
	match($0, /([0-9]+)$/, total_out);
	next;
    }

    if ($0 ~ /INFO regex matches read1:/)
    {
	match($0, /([0-9]+)$/, raw_match);
	next;
    }

    if ($0 ~ /INFO regex does not match read1:/)
    {
	match($0, /([0-9]+)$/, raw_unmatch);
	next;
    }

    if ($0 ~ /INFO Filtered cell barcode. Not correctable:/)
    {
	match($0, /([0-9]+)$/, filtered);
	next;
    }

    if ($0 ~ /INFO False cell barcode. Error-corrected:/)
    {
	match($0, /([0-9]+)$/, corrected);
	next;
    }
}
END {
    print sample[1], total_in[1], total_out[1], raw_match[1], raw_unmatch[1], filtered[1], corrected[1];
}
