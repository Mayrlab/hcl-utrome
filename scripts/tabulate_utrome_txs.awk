#!/usr/bin/awk -f
BEGIN {
    FS="\t";
    OFS=",";
    gencode=0;
    upstream=0;
    downstream=0;
} {
    if ($3 ~ /transcript/) {
	match($9, /ID "[^-"]+(-UTR)?(-|\+)?[0-9]*"/, id);
	switch (id[2]) {
	    case "+":
		downstream++;
		break;
	    case "-":
		upstream++;
		break;
	    default:
		gencode++;
	}
    }
} END {
    print "site_type", "count";
    print "gencode", gencode;
    print "upstream", upstream;
    print "downstream", downstream;
}
