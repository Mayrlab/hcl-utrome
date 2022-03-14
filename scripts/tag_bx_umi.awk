#!/usr/bin/awk -f

BEGIN {
    OFS="\t";
}
{
    if ($1 ~ /^@/) {
	print;
    } else {
	split($1, ids, "_");
	$1="";
	print ids[1]$0,"CB:Z:"ids[2],"RX:Z:"ids[3];
    }
}
