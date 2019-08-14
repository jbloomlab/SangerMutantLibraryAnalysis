"""Script to analyze mutations to a coding sequence.

Written by Jesse Bloom, 2013

Edited by Hugh Haddox, October-16-2015
Added optional command line argument parsing by Mike Doud October 23 2015
Fixed x-axis integer bug...
Edited by Allie Greaney, November-22-2018
Fixed x-axis bar alignment bug and converted to python3.
Edited by Allie Greaney, August-13-2019
Fixed bug that led to non-counting of STOP codons.
"""

import re
import os
import time
import math
import random
import matplotlib
matplotlib.use('pdf') # use the PDF backend
matplotlib.rc('legend', fontsize=12)
import pylab
import scipy.stats
import argparse





def TranslateCodon(codon):
    """Returns one-letter amino acid code for *codon*.

    *codon* is a 3-letter string giving a valid codon."""
    genetic_code = {'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L',
        'CTA':'L', 'CTG':'L', 'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M', 'GTT':'V',
        'GTC':'V', 'GTA':'V', 'GTG':'V', 'TCT':'S', 'TCC':'S', 'TCA':'S',
        'TCG':'S', 'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'ACT':'T',
        'ACC':'T', 'ACA':'T', 'ACG':'T', 'GCT':'A', 'GCC':'A', 'GCA':'A',
        'GCG':'A', 'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*',
        'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q', 'AAT':'N', 'AAC':'N',
        'AAA':'K', 'AAG':'K', 'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
        'TGT':'C', 'TGC':'C', 'TGA':'*', 'TGG':'W', 'CGT':'R',
        'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGT':'S', 'AGC':'S', 'AGA':'R',
        'AGG':'R', 'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}
    return genetic_code[codon.upper()]


def PlotMutationClustering(mutation_nums_by_clone, ncodons, plotfile, title, mutstart, nsimulations=1000):
    """Plots clustering of mutations versus null expectation of no clustering.

    This function addresses the question of whether clones with multiple
    mutations tend to have those mutations clustered in primary sequence.
    Clones with <2 mutations are not considered. For every clone with >= 2
    mutations, records the distance in primary sequence between each pair
    of mutations. Then performs nsimulations simulations of placing this
    number of mutations at random on a gene, and records the distance in
    primary sequence between each pair in the simulated sequences. Plots
    the actual distributions of primary sequence distances versus
    the simulated distribution.

    mutation_nums_by_clone -> list, with an entry for each clone. Each of
        these entries is itself a list. The entries in these sublists are
        numbers (1, 2, 3, ...) of the codon positions mutated in that clone.
    ncodons -> integer number of codons in the gene.
    plotfile -> name of the plot file we create.
    title -> string giving the plot title.
    nsimulations -> number of simulations for each clone. Is 1000 by default.
    mutstart -> specifies the first codon in the mutated segment of the gene
        (integer). This variable will be used to truncate the gene to the
        appropriate length for simulating the random distribution of distances
        between mutations.
    """
    actual_distances = dict([(i, 0) for i in range(1, ncodons - mutstart + 1)])
    simulated_distances = dict([(i, 0) for i in range(1, ncodons - mutstart + 1)])
    codons = [i for i in range(1, ncodons - mutstart + 2)]
    nactual = nsimulated = 0
    for mutpositions in mutation_nums_by_clone:
        nmuts = len(mutpositions)
        if nmuts < 2:
            continue
        for i in range(nmuts):
            for j in range(i + 1, nmuts):
                d = abs(mutpositions[i] - mutpositions[j])
                actual_distances[d] += 1
                nactual += 1
        for isimulate in range(nsimulations):
            simulpositions = random.sample(codons, nmuts)
            for i in range(nmuts):
                for j in range(i + 1, nmuts):
                    d = abs(simulpositions[i] - simulpositions[j])
                    simulated_distances[d] += 1
                    nsimulated += 1
    actual_cumul = []
    simulated_cumul = []
    actual_tot = simul_tot = 0.0
    for d in range(1, ncodons - mutstart + 1):
        actual_tot += actual_distances[d] / float(nactual)
        simul_tot += simulated_distances[d] / float(nsimulated)
        actual_cumul.append(actual_tot)
        simulated_cumul.append(simul_tot)
    pylab.figure(figsize=(4.5, 2.25))
    (lmargin, rmargin, bmargin, tmargin) = (0.13, 0.01, 0.21, 0.07)
    pylab.axes([lmargin, bmargin, 1.0 - lmargin - rmargin, 1.0 - bmargin - tmargin])
    barwidth = 0.7
    xs = [x for x in range(1, ncodons - mutstart + 1)]
    assert len(xs) == len(actual_cumul) == len(simulated_cumul)
    pred = pylab.plot(xs, simulated_cumul, 'b--')
    actual = pylab.plot(xs, actual_cumul, 'r-')
    pylab.gca().set_xlim([0, ncodons])
    pylab.gca().set_ylim([0, 1])
    pylab.gca().xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(6))
    pylab.gca().yaxis.set_major_locator(matplotlib.ticker.FixedLocator([0, 0.5, 1]))
    pylab.xlabel('distance between pairs of mutations')
    pylab.ylabel('cumulative fraction')
    pylab.legend((actual[0], pred[0]), ('actual', 'expected'), loc='lower right', numpoints=1, handlelength=2, ncol=2, borderaxespad=0.4, handletextpad=0.4, columnspacing=1.1)
    pylab.title(title, fontsize=12)
    pylab.savefig(plotfile)
    time.sleep(0.5)
    pylab.show()


def PlotCodonMutNTComposition(allmutations, plotfile, title):
    """Plots nucleotide composition of mutant codons and mtuated.

    For each site containing a mutated codon, looks at the nucleotide
    composition of all three sites at the original and new codon. Plots the
    overall frequency of each nucleotide at these sites.

    allmutations -> list of all mutations as tuples (wtcodon, r, mutcodon)
    plotfile -> name of the plot file we create.
    title -> string giving the plot title.
    """
    nts = ['A', 'T', 'C', 'G']
    wtntcounts = dict([(nt, 0) for nt in nts])
    mutntcounts = dict([(nt, 0) for nt in nts])
    ntot = float(3 * len(allmutations))
    for (wtcodon, r, mutcodon) in allmutations:
        for nt in wtcodon:
            wtntcounts[nt] += 1 / ntot
        for nt in mutcodon:
            mutntcounts[nt] += 1 / ntot
    pylab.figure(figsize=(3.5, 2.25))
    (lmargin, rmargin, bmargin, tmargin) = (0.16, 0.01, 0.21, 0.07)
    pylab.axes([lmargin, bmargin, 1.0 - lmargin - rmargin, 1.0 - bmargin - tmargin])
    barwidth = 0.35
    pylab.gca().set_xlim([-0.5, 3.5])
    pylab.gca().set_ylim([0, max(list(wtntcounts.values()) + list(mutntcounts.values())) * 1.35])
    xs = [i for i in range(len(nts))]
    nwt = pylab.bar([x - barwidth/2 for x in xs], [wtntcounts[nt] for nt in nts], width=barwidth, color='blue')
    nmut = pylab.bar([x + barwidth/2 for x in xs], [mutntcounts[nt] for nt in nts], width=barwidth, color='red')
    #pred = pylab.plot(xs, nexpected, 'rx', markersize=6, mew=3)
    #pylab.gca().xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
    pylab.gca().yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5))
    pylab.xlabel('nucleotide')
    pylab.ylabel('codon composition')
    pylab.legend((nwt[0], nmut[0]), ('parent', 'mutant'), loc='upper center', numpoints=1, handlelength=1, ncol=2, borderaxespad=0.01, columnspacing=1.1, handletextpad=0.4)
    pylab.title(title, fontsize=12)
    pylab.xticks(pylab.arange(0, 4, 1), tuple(nts))
    pylab.savefig(plotfile)
    time.sleep(0.5)
    pylab.show()


def PlotNCodonMuts(allmutations, plotfile, title):
    """Plots number of nucleotide changes per codon mutation.

    allmutations -> list of all mutations as tuples (wtcodon, r, mutcodon)
    plotfile -> name of the plot file we create.
    title -> string giving the plot title.
    """
    pylab.figure(figsize=(3.5, 2.25))
    (lmargin, rmargin, bmargin, tmargin) = (0.16, 0.01, 0.21, 0.07)
    pylab.axes([lmargin, bmargin, 1.0 - lmargin - rmargin, 1.0 - bmargin - tmargin])
    nchanges = {1:0, 2:0, 3:0}
    nmuts = len(allmutations)
    for (wtcodon, r, mutcodon) in allmutations:
        assert 3 == len(wtcodon) == len(mutcodon)
        diffs = len([i for i in range(3) if wtcodon[i] != mutcodon[i]])
        nchanges[diffs] += 1
    barwidth = 0.6
    xs = [1, 2, 3]
    nactual = [nchanges[x] for x in xs]
    nexpected = [nmuts * 9. / 63., nmuts * 27. / 63., nmuts * 27. / 63.]
    bar = pylab.bar([x for x in xs], nactual, width=barwidth)
    pred = pylab.plot(xs, nexpected, 'rx', markersize=6, mew=3)
    pylab.gca().set_xlim([0.5, 3.5])
    pylab.gca().set_ylim([0, max(nactual + nexpected) * 1.1])
    pylab.gca().xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    pylab.gca().yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5, integer=True))
    pylab.xlabel('nucleotide changes in codon')
    pylab.ylabel('number of mutations')
    pylab.legend((bar[0], pred[0]), ('actual', 'expected'), loc='upper left', numpoints=1, handlelength=0.9, borderaxespad=0, handletextpad=0.4)
    pylab.title(title, fontsize=12)
    pylab.savefig(plotfile)
    time.sleep(0.5)
    pylab.show()




def PlotNMutDist(nmutations, plotfile, title):
    """Plots number of mutations per gene versus Poisson distribution.

    nmutations -> list of number of mutations per gene.
    plotfile -> name of the plot file that we create.
    title -> string giving the plot title.
    """
    pylab.figure(figsize=(3.5, 2.25))
    (lmargin, rmargin, bmargin, tmargin) = (0.20, 0.01, 0.21, 0.07)
    pylab.axes([lmargin, bmargin, 1.0 - lmargin - rmargin, 1.0 - bmargin - tmargin])
    nseqs = len(nmutations)
    mavg = scipy.mean(nmutations)
    barwidth = 0.8
    xmax = max(nmutations) + 2
    nmuts = [i for i in range(xmax + 1)]
    nactual = [nmutations.count(n) for n in nmuts]
    npoisson = [nseqs * math.exp(-mavg) * mavg**n / math.factorial(n) for n in nmuts]
    bar = pylab.bar([n for n in nmuts], nactual, width=barwidth)
    pred = pylab.plot(nmuts, npoisson, 'rx', markersize=6, mew=3)
    pylab.gca().set_xlim([-0.5, xmax + 0.5])
    pylab.gca().set_ylim([0, max(nactual) + 1.5])
    pylab.gca().xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    pylab.gca().yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    pylab.xlabel('number of mutated codons')
    pylab.ylabel('number of clones')
    pylab.legend((bar[0], pred[0]), ('actual', 'Poisson'), loc='upper right', numpoints=1, handlelength=1.2, ncol=1, borderaxespad=0)
    #
    # Kolmogorov-Smirnov test
    def F(n):
        return scipy.stats.poisson.cdf(n, mavg)
    (d, p) = scipy.stats.kstest(scipy.array(nmutations), F)
#    pylab.text(0.98, 0.6,
#        'P = %.3f for Kolmogorov-Smirnov test of\n' % p +\
#        'whether actual distribution would differ\n' +\
#        'from Poisson at least this much by chance.',
#        horizontalalignment='right',
#        verticalalignment='bottom',
#        transform=pylab.gca().transAxes,
#        fontsize=11)
    pylab.title(title, fontsize=12)
    pylab.savefig(plotfile)
    time.sleep(0.5)
    pylab.show()


def PlotGeneMutDist(genelength, sub_nums, indel_nums, plotfile, cumulplotfile, title, mutstart):
    """Plots mutation distribution along a gene.

    genelength -> length of gene (integer)
    sub_nums -> list of positions of substitutions (1, 2, ... numbering)
    indel_nums -> list of positions of indels (1, 2, ... numbering)
    plotfile -> string giving plot file name for plot with lines at each
        site.
    cumulplotfile -> string giving name of cumulative distribution plot file.
    title -> string giving plot title
    mutstart -> specifies the first codon in the mutated portion of the gene (integer).
        The cumulative and uniform distributions will begin at this position.
    """
    nsubs = len(sub_nums)
    if not nsubs:
        raise ValueError("empty sub_nums")
    xs = [i for i in range(mutstart, genelength + 1)]
    subs = dict([(x, 0) for x in xs])
    indels = dict([(x, 0) for x in xs])
    for x in sub_nums:
        subs[x] += 1
    for x in indel_nums:
        indels[x] += 1
    subs = [subs[x] for x in xs]
    indels = [-indels[x] for x in xs]
    barwidth = 1.0
    xlefts = [x - barwidth / 2. for x in xs]
    pylab.figure(figsize=(5.5, 2.5))
    (lmargin, rmargin, bmargin, tmargin) = (0.1, 0.03, 0.2, 0.09)
    pylab.axes([lmargin, bmargin, 1.0 - lmargin - rmargin, 1.0 - bmargin - tmargin])
    pylab.bar(xlefts, subs, width=barwidth)
    pylab.bar(xlefts, indels, width=barwidth)
    pylab.gca().set_xlim([mutstart - 0.5, genelength + 0.5])
    yticker = matplotlib.ticker.MultipleLocator(1)
    pylab.gca().yaxis.set_major_locator(yticker)
    ymax = max(subs + indels)
    pylab.gca().set_ylim([-ymax, ymax])
    pylab.title(title, fontsize=12)
    pylab.xlabel('codon number')
    pylab.ylabel('indels <-> subst.')
    pylab.savefig(plotfile)
    time.sleep(0.5)
    pylab.show()
    # make cumulative distribution plot
    cumul = []
    for i in range(len(subs)):
        if cumul:
            cumul.append(cumul[-1] + subs[i] / float(nsubs))
        else:
            cumul.append(subs[i] / float(nsubs))
    pylab.figure(figsize=(4.5, 2.25))
    (lmargin, rmargin, bmargin, tmargin) = (0.13, 0.04, 0.21, 0.07)
    pylab.axes([lmargin, bmargin, 1.0 - lmargin - rmargin, 1.0 - bmargin - tmargin])
    linear = [(x-mutstart+1) / float(genelength-mutstart+1) for x in xs]
    pylab.plot(xs, cumul, 'r-', label='actual')
    pylab.plot(xs, linear, 'b--', label='uniform')
    pylab.gca().set_xlim([mutstart, genelength])
    pylab.gca().set_ylim([0, 1])
    pylab.ylabel('cumulative fraction')
    pylab.gca().yaxis.set_major_locator(matplotlib.ticker.FixedLocator([0, 0.5, 1]))
    pylab.xlabel('codon number')
    pylab.legend(loc='upper left')
    pylab.title(title, fontsize=12)
    pylab.savefig(cumulplotfile)
    time.sleep(0.5)
    pylab.show()
    pylab.clf()
    pylab.close()


def ReadClone(line, seq):
    """Parses information from a line specifying mutations for a clone.

    line -> line containing mutation information
    seq -> the wildtype sequence.
    Returns the following tuple: (name, mutations, indels)
        name -> string giving name of clone
        mutations -> list of 3-tuples of wildtype codon, codon number, mutant codon
        indels -> list of codon numbers of sites where indels occur
    """
    submatch = re.compile('^(?P<wt>[ATCG]+)(?P<num>\d+)(?P<mut>[ATCG]+)$')
    indelmatch = re.compile('(ins|del)([ATCG]+)(?P<num>\d+)')
    (name, muts) = line.split(':')
    muts = muts.split(',')
    if len(muts) == 1 and muts[0].strip() == 'None':
        return (name, [], [])
    mutations = []
    indels = []
    for mutstring in muts:
        mutstring = mutstring.strip()
        if not mutstring:
            raise ValueError("No mutations specified for a clone. If a clone has no mutations, enter None in the mutation list.")
        m = submatch.search(mutstring)
        if m:
            (wt, num, mut) = (m.group('wt'), int(m.group('num')), m.group('mut'))
            assert len(wt) == len(mut)
            if seq[num - 1 : num - 1 + len(wt)] != wt:
                raise ValueError("Mismatch at %s. Make sure you entered this mutation correctly -- probably you entered the wrong wildtype identity or misnumbered the mutation." % mutstring)
            (icodon, ipos) = divmod(num - 1, 3)
            icodon = icodon + 1
            if ipos + len(wt) > 3:
                raise ValueError("Mutation spans multiple codons")
            wtcodon = seq[3 * (icodon - 1) : 3 * icodon]
            assert wt in wtcodon
            mutcodon = list(wtcodon)
            for i in range(len(mut)):
                mutcodon[i + ipos] = mut[i]
            mutations.append((wtcodon, icodon, mutcodon))
        else:
            m = indelmatch.search(mutstring)
            if not m:
                raise ValueError("Couldn't match sub or indel: %s" % mutstring)
            num = int(m.group('num'))
            (icodon, ipos) = divmod(num - 1, 3)
            icodon = icodon + 1
            indels.append(icodon)
    return (name, mutations, indels)


def main():
    """Main body of script."""

    parser = argparse.ArgumentParser()
    parser.add_argument("--outputprefix", help="optional prefix for output files generated")
    parser.add_argument("--seqfile", help="name of FASTA file containing the gene sequence")
    parser.add_argument("--mfile", help="name of the file containing the list of mutations")
    parser.add_argument("--mutstart", help="position of the first codon in the mutated segment of the gene")
    parser.add_argument("--title", help="title for plots generated")
    args = parser.parse_args()

    print("\nBeginning analysis.")

    if not args.seqfile:
        seqfile = input("\nEnter the name of the FASTA file containing the gene sequence: ").strip()
    else:
        seqfile = args.seqfile
    if not os.path.isfile(seqfile):
        raise IOError("Cannot find specified file %s" % seqfile)
    seq = open(seqfile).readlines()[1].strip().upper()
    print("Read a coding sequence of length %d" % len(seq))
    assert len(seq) % 3 == 0, "Sequence length is not a multiple of three."
    ncodons = len(seq) // 3

    # read sequence
    if not args.mfile:
        mfile = input("\nEnter the name of the file containing the list of mutations: ").strip()
    else:
        mfile = args.mfile
    if not os.path.isfile(mfile):
        raise IOError("Cannot find specified file %s" % mfile)

    # get the position of the first codon in the mutated portion of the gene
    if not args.mutstart:
        mutstart = int(input("\nEnter the position of the first codon in the mutated segment of the gene: ").strip())
    else:
        mutstart = int(args.mutstart)

    # begin looping over input libraries
    print("\nReading mutations from %s" % mfile)
    clones = [line for line in open(mfile).readlines() if (not line.isspace()) and line[0] != '#']
    print("Read entries for %d clones" % len(clones))
    clone_d = {}
    sub_nums = []
    indel_nums = []
    nmutations = []
    allmutations = []
    mutation_nums_by_clone = []
    for line in clones:
        (name, mutations, indels) = ReadClone(line, seq)
        if name in clone_d:
            raise ValueError("duplicate clone %s" % name)
        clone_d[name] = True
        nmutations.append(len(mutations))
        allmutations += mutations
        imutation_nums_by_clone = []
        for (wtcodon, icodon, mutcodon) in mutations:
            if icodon < mutstart:
                raise ValueError("This line reports a mutation before the beginning of the mutated segment of the gene: %s" %line)
            imutation_nums_by_clone.append(icodon)
            sub_nums.append(icodon)
        mutation_nums_by_clone.append(imutation_nums_by_clone)
        for icodon in indels:
            if icodon < mutstart:
                raise ValueError("This line reports an indel before the beginning of the mutated segment of the gene: %s" %line)
            indel_nums.append(icodon)
    sub_nums.sort()
    indel_nums.sort()
    print("\nSubstitutions begin at following positions: %s" % ', '.join([str(s) for s in sub_nums]))
    print("\nIndels begin at following positions: %s" % ', '.join([str(s) for s in indel_nums]))
    denom = float((ncodons - mutstart + 1) * len(clones))
    print("\nFound a total of %d substitutions out of %d codons sequenced (%.4f)" % (len(allmutations), (ncodons - mutstart + 1) * len(clones), len(allmutations) / denom))
    n_nmuts = {1:0, 2:0, 3:0}
    n_muttypes = {'synonymous':0, 'nonsynonymous':0, 'stop codon':0}
    for (wt, i, m) in allmutations:
        ndiffs = len([i for i in range(3) if wt[i] != m[i]])
        n_nmuts[ndiffs] += 1
        m = ''.join(m)
        wtaa = TranslateCodon(wt)
        mutaa = TranslateCodon(m)
        if wtaa == mutaa:
            n_muttypes['synonymous'] += 1
        elif wtaa != '*' and mutaa == '*':
            n_muttypes['stop codon'] += 1
        else:
            n_muttypes['nonsynonymous'] += 1
    print("\nHere are the fractions with different numbers of nucleotide mutations:")
    for n in range(1, 4):
        print("  %d nucleotide mutations: %.5f" % (n, n_nmuts[n] / denom))
    print("\nHere are the fractions of mutation types")
    for key in ['synonymous', 'nonsynonymous', 'stop codon']:
        print("  %s: %.5f" % (key, n_muttypes[key] / denom))
    nclones = len(clone_d)
    print("\nOverall summary:\n%d clones, avg. %.1f codon substitutions, avg. %.1f indels" % (nclones, len(sub_nums) / float(nclones), len(indel_nums) / float(nclones)))
    print("\nNow creating the output PDF plot files...")

    if not args.title:
        title = ''
    else:
        title = args.title

    if not args.outputprefix:
        outputprefix = ''
    else:
        outputprefix = args.outputprefix

    PlotGeneMutDist(ncodons, sub_nums, indel_nums, "%s_mutpositions.pdf" % outputprefix, "%s_mutpositions_cumulative.pdf" % outputprefix, title, mutstart)
    os.system('convert -density 150 %s_mutpositions.pdf %s_mutpositions.jpg' % (outputprefix, outputprefix))
    os.system('convert -density 150 %s_mutpositions_cumulative.pdf %s_mutpositions_cumulative.jpg' % (outputprefix, outputprefix))
    PlotNMutDist(nmutations, '%s_nmutdist.pdf' % outputprefix, '')
    os.system('convert -density 150 %s_nmutdist.pdf %s_nmutdist.jpg' % (outputprefix, outputprefix))
    PlotNCodonMuts(allmutations, '%s_ncodonmuts.pdf' % outputprefix, '')
    os.system('convert -density 150 %s_ncodonmuts.pdf %s_ncodonmuts.jpg' % (outputprefix, outputprefix))
    PlotCodonMutNTComposition(allmutations, '%s_codonmutntcomposition.pdf' % outputprefix, '')
    os.system('convert -density 150 %s_codonmutntcomposition.pdf %s_codonmutntcomposition.jpg' % (outputprefix, outputprefix))
    PlotMutationClustering(mutation_nums_by_clone, ncodons, '%s_mutationclustering.pdf' % outputprefix, '', mutstart)
    os.system('convert -density 150 %s_mutationclustering.pdf %s_mutationclustering.jpg' % (outputprefix, outputprefix))
    print("The output PDF file plots have now all been created.\n\nScript complete.")


main() # run the script
