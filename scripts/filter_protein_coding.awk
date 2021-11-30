#!/usr/bin/env awk -f
BEGIN {
    FS="\t";
}
{
    if ($3 ~ /gene/) {
	## print pc genes
	if ($9 ~ /gene_type "protein_coding"/) {
	    print
	}
    } else {
	## print all pc transcripts (and children)
	if ($9 ~ /transcript_type "protein_coding"/) {
	    print
	}
    }

    ## retain header
    if ($1 ~ /^#/) { print }
}
