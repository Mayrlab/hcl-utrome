#!/usr/bin/env awk -f

BEGIN { OFS=","; }
{
    if ($0 ~ /Total reads:/)
    {
        match($0, /([0-9]+)$/, total);
        next;
    }

    if ($0 ~ /Aligned 0 time:/)
    {
	match($0, /: ([0-9]+) \(([^%]+)/, unaligned);
	next;
    }
    
    if ($0 ~ /Aligned 1 time:/)
    {
	match($0, /: ([0-9]+) \(([^%]+)/, unique);
	next;
    }
    
    if ($0 ~ /Aligned >1 times:/)
    {
	match($0, /: ([0-9]+) \(([^%]+)/, multi);
	next;
    }
}
END {
    match(FILENAME, /([^/]+)\.assembled\.log$/, sample);
    print sample[1], total[0], unaligned[1], unique[1], multi[1];
}
