#!/usr/bin/python



import pysam

import argparse

from natsort import humansorted

from collections import defaultdict

import numpy as np





def parse_args():

    parser = argparse.ArgumentParser(description='Parse BAM component filter arguments.')

    parser.add_argument('--bam_file', '-b', type=str, required=True,

                        help='Path to the BAM file(s) to parse. (comma-separated if multiple)')

    parser.add_argument('--out_file', '-o', required=True, type=str, default='counts.tsv',

                        help='Pat of file to write, a TSV of format <contig> <contig> <link_count>. Default: counts.tsv')

    parser.add_argument('--out_fmt', '-f', required=False, type=str, default="counts",

                        help='Format of file to write. Default: counts (name1 name2 counts). others: index_counts (index1 index2 counts)' \

                             'npy (NYI), matrix (NYI)')

    args = parser.parse_args()

    return vars(args)





# parse bam file to links

def parse_bam_to_link_counts(bamfiles):

    '''Parse the bam file to make it into a Hi-C graph among contig nodes with read edges. (Edges are unweighted).

    Also pull length info for the contigs to filter on later.



    Args:

        bamfiles ([str]): List of bamfile paths to parse (can be length 1).



    Returns:

        {str: {str: bool}}: Representation of the Hi-C graph.

        {str: int}: Mapping of contig names to their sequence length.



    '''

    net = defaultdict(lambda: defaultdict(int))

    if len(bamfiles) == 0:

        RuntimeError("No BAM files provided!")

    print "parsing bam(s)"

    for bamfile in bamfiles:

        bam = pysam.AlignmentFile(bamfile, 'rb')

        num = 0

        num_dupes = 0

        is_first_read = True

        refs = bam.references

        lengths = bam.lengths

        ref_lens = {}

        ref_reads_counts = {}

    mate = ''

    for read in bam:

        if is_first_read:

            num += 1

            is_first_read = False

            read1_id = read.query_name

            if read.is_duplicate:

                num_dupes += 1

                continue

            mate = [read.next_reference_start, read.next_reference_id]

            ref1 = refs[read.reference_id]

            ref_lens[ref1] = lengths[read.reference_id]

            ref2 = refs[read.next_reference_id]

            ref_lens[ref2] = lengths[read.next_reference_id]

            if num % 10000000 == 0:

                print "parsed {0} read pairs from file {1}, {2} filtered out as duplicates".format(num, bamfile,

                                                                                                   num_dupes)

        else:

            are_paired = read.query_name == read1_id

            if are_paired:

                net[ref1][ref2] += 1

                net[ref2][ref1] += 1

                if ref1 in ref_reads_counts:

                    ref_reads_counts[ref1] += 1

                else:

                    ref_reads_counts[ref1] = 0

                if ref2 in ref_reads_counts:

                    ref_reads_counts[ref2] += 1

                else:

                    ref_reads_counts[ref2] = 0

                is_first_read = True

            else:

                is_first_read = False

                # print read.query_name, read1_id

                read1_id = read.query_name

                if read.is_duplicate:

                    num_dupes += 1

                    continue

                mate = [read.next_reference_start, read.next_reference_id]

                ref1 = refs[read.reference_id]

                ref_lens[ref1] = lengths[read.reference_id]

                ref2 = refs[read.next_reference_id]

                ref_lens[ref2] = lengths[read.next_reference_id]

                if num % 10000000 == 0:

                    print "parsed {0} read pairs from file {1}, {2} filtered out as duplicates".format(num, bamfile,

                                                                                                       num_dupes)

            # print "M02014:69:000000000-BT93H:1:1101:1721:15359" is "M02014:69:000000000-BT93H:1:1101:1721:15359"



            # raise ValueError("BAM file is not sorted by read name!! (samtools sort -n)")

            if read.query_name != read1_id or mate != [read.reference_start, read.reference_id]:

                RuntimeError(

                    "BAM file is not sorted by read name!! (samtools sort -n). Alternately, filtering may have removed reads rather than read pairs?")

    print "parsed {0} read pairs from file {1}, filtered out {2} duplicates".format(num, bamfile, num_dupes)



    if num_dupes == 0:

        print "number of duplicates is zero- this may be because the duplicates flag was not set in BAM (e.g. by samblaster)"



    return net, ref_lens, ref_reads_counts





def write_reads_len_on_contig(ref_reads_counts, ref_lens, outfile):

    '''Write out all those number of mapped reads on contigsto a file.'''

    with open(outfile+"_stat_contigs", 'w') as out:

        for contig in humansorted(ref_reads_counts.keys()):

            out.write(contig +"\t" +str(ref_lens[contig])+"\t" +str(ref_reads_counts[contig]) + "\n")





def write_net_as_counts(net, outfile):

    '''Write out all those counts from the net dict of defaultdicts to a file.'''

    print "writing counts data to {0}".format(outfile)
    s = set()

    with open(outfile, 'w') as out:

        for contig1 in humansorted(net.keys()):

            this_dict = net[contig1]

            for contig2 in humansorted(this_dict.keys()):

                r1 = contig1+contig2
                r2 = contig2+contig1
                if(r1 not in s and r2 not in s):
                    s.add(r1)
                    s.add(r2)
                    # print contig1, contig2, this_dict[contig2]
                    outstr = "\t".join([contig1, contig2, str(this_dict[contig2])])
                    outstr = "\t".join([contig2, contig1, str(this_dict[contig2])])
                    out.write(outstr + "\n")





def write_net_as_counts_index(net, outfile):

    '''Write out all those counts from the net dict of defaultdicts to a file.

    Using index names to accommodate tools like EVR that expect such.'''

    print "writing counts data to {0}".format(outfile)

    keys = humansorted(net.keys())

    with open(outfile, 'w') as out:

        for contig1 in keys:

            idx_print1 = str(keys.index(contig1))

            print contig1

            this_dict = net[contig1]

            for contig2 in keys:

                idx_print2 = str(keys.index(contig2))

                # print contig1, contig2, this_dict[contig2]

                outstr = "\t".join([idx_print1, idx_print2, str(this_dict[contig2])])

                out.write(outstr + "\n")

    # write out contig names separately for tools that require indices...

    with open(outfile + ".names", "w") as outnames:

        outnames.write("\n".join(keys) + "\n")





def write_net_as_npy(net, outfile):

    '''Write out all those counts from the net dict as an npy matrix that can be used by e.g. pastis.

    Requires taking numpy dependency, may or mayn't be worth it.

    '''

    UserWarning("zeroing out diagonal!!!! pastis requires this")

    contigs = humansorted(net.keys())

    counts = np.zeros(shape=(len(contigs), len(contigs)), dtype=int, order="F")

    keys = humansorted(net.keys())

    for contig1 in keys:

        # print contig1

        idx1 = keys.index(contig1)

        this_dict = net[contig1]

        for contig2 in this_dict.keys():

            if contig1 == contig2:

                continue

            idx2 = keys.index(contig2)

            # print contig1, contig2, this_dict[contig2]

            counts[idx1, idx2] = this_dict[contig2]

    counts.dump(outfile)

    np.save(outfile, counts)





def write_net_as_matrix(net, outfile):

    '''Write out all those counts from the net dict as a flat matrix file'''

    pass





def read_net_from_counts(counts_file):

    '''Read a previously parsed network from a counts file. Assumes a 3-column file. Not used when run as a script.'''

    print "reading counts into a net from the file {0}.".format(counts_file)

    with open(counts_file) as file:

        net = defaultdict(lambda: defaultdict(int))

        for line in file:

            fields = line.strip().split()

            ref1 = fields[0]

            ref2 = fields[1]

            count = fields[2]

            net[ref1][ref2] = count

            net[ref2][ref1] = count

    return net





def main():

    c_args = parse_args()

    bams = c_args["bam_file"].split(",")

    net, ref_lens, ref_reads_counts = parse_bam_to_link_counts(bams)

    print "writing out data"

    if c_args["out_fmt"] == "counts":

        write_net_as_counts(net=net, outfile=c_args["out_file"])

    elif c_args["out_fmt"] == "index_counts":

        write_net_as_counts_index(net=net, outfile=c_args["out_file"])

    elif c_args["out_fmt"] == "npy":

        write_net_as_npy(net=net, outfile=c_args["out_file"])

    write_reads_len_on_contig(ref_reads_counts=ref_reads_counts, ref_lens = ref_lens, outfile=c_args["out_file"])





if __name__ == "__main__":

    main()


