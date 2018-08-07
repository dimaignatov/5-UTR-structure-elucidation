import re
import os
import HTSeq
import pickle
import numpy
from Bio import SeqIO
from Bio.Seq import Seq

def write_pickle(ga, pickle_fn):
    ''' Generates pickle file for a Genomic Array object.

Requires the Genomic Array object and the name for pickle file
'''
    cw = open(pickle_fn, 'w')
    pickle.dump(ga, cw)
    cw.close()

def write_artemis(ga, artemis_fn):
    ''' Generates .artemis file with profile.

Requires Genomic array object with coverage or mismatch data, and the filename
for .artemis file.
'''
    genome_id = ga.chrom_vectors.keys()[0]
    plus_gi = ga.chrom_vectors[genome_id]['+'].iv
    minus_gi = ga.chrom_vectors[genome_id]['-'].iv
    plus = list(ga[plus_gi])
    minus = list(ga[minus_gi])

    fw = open(artemis_fn, 'w')
    for (P, M) in zip(plus, minus):
        fw.write('{0}\t{1}\n'.format(P, M))
    fw.close()

def write_bed(ga, bed_plus_fn, bed_minus_fn):
    ''' Generates .bedgraph files for plus and minus genome strands.

Requires Genomic array object with coverage or mismatch data, and the filenames
for these .bedgraph files.
'''
    genome_id = ga.chrom_vectors.keys()[0]
    id_mod = re.sub(r'\.\d+', '', genome_id)

    ga.write_bedgraph_file('temp.bedgraph', '+', track_options="color=255,0,0")
    BED_PLUS = open(bed_plus_fn, 'w')
    with open('temp.bedgraph') as TEMP:
        for line in TEMP:
            if not genome_id in line:
                BED_PLUS.write(line)
                continue
            line = line.rstrip()
            genome_id, start, end, val = re.split(r'\t', line)
            BED_PLUS.write('{0}\t{1}\t{2}\t{3}\n'.format(id_mod, start, end, float(val)))
    BED_PLUS.close()
    os.remove('temp.bedgraph')


    ga.write_bedgraph_file('temp.bedgraph', '-')
    BED_MINUS = open(bed_minus_fn, 'w')
    with open('temp.bedgraph') as TEMP:
        for line in TEMP:
            if not genome_id in line:
                BED_MINUS.write(line)
                continue
            line = line.rstrip()
            genome_id, start, end, val = re.split(r'\t', line)
            BED_MINUS.write('{0}\t{1}\t{2}\t{3}\n'.format(id_mod, start, end, -float(val)))
    BED_MINUS.close()
    os.remove('temp.bedgraph')

def parse_md(md):
    ''' Parses MD field in SAM file.

MD field should be in a raw format, e.g. '25G45^A0C28'.
Returns a vector of tuples containing a relative 0-based coordinate of
the substitution and the original base in the genome. Does not search
for indels. In case of none substitutions, returns an empty vector.
'''

    substitutions = []
    position = 0
    p1 = re.compile(r'\d+\D*')
    p2 = re.compile(r'(\d+)(\D*)')
    p3 = re.compile(r'\^')
    fields = p1.findall(md)
    for field in fields:
        m = p2.match(field)
        offset = m.group(1)
        position = position + int(offset)
        base = m.group(2)
        if not base:
            return substitutions
        if p3.match(base):
            position = position + len(base) - 1
            next
        else:
            substitutions.append((position, base))
            position = position + 1

def generateProfiles(sam_fn, fasta_fn='NC_003210.1.fa'):
    ''' Creates coverage and mismatch profiles for DMS-MaPseq sample.

Requires the .sam file of Bowtie2 alignment of DMS-MaPseq reads to the bacterial genome.
Also requires the .fasta file of the genome, to which the reads were aligned.
Creates the profiles for coverage and number of mismatches for each nucleotide
in the genome. The profiles are saved as .pickle files in the form of GenomicArray
objects (See HTSeq library docs). Additionally the profiles can be visualized with
Artemis and IGV genome browsers.
'''
    # read the genome sequence from fasta file
    genome = SeqIO.read(fasta_fn, "fasta")
    genome_length = len(genome.seq)
    # create genomic arrays to store coverage and mismatch data
    cvg = HTSeq.GenomicArray({genome.id: genome_length}, stranded=True,
                             typecode="i")
    mis = HTSeq.GenomicArray({genome.id: genome_length}, stranded=True,
                             storage='ndarray', typecode="i")
    # read the paired-end sam file
    sam_reader = HTSeq.SAM_Reader(sam_fn)
    i = 0
    for first, second in HTSeq.pair_SAM_alignments(sam_reader):
        i += 1
        if not i % 100000:
            print sam_fn, '->', i
        if not (first.proper_pair and first.proper_pair):
            continue
        
        # Add fragment coverage to the coverage profile 
        second.iv.strand = first.iv.strand # The first read determines the fragment strand
        # If first and second reads overlap, the coverage is calculated for the whole fragment
        if first.iv.overlaps(second.iv):
            first.iv.extend_to_include(second.iv)
            cvg[first.iv] += 1
        else: # Alternatively, the coverage is calculated for each read separately
            cvg[first.iv] += 1
            cvg[second.iv] += 1

        # Add unique mismatches from every pair of reads to the mismatch profile
        mism_1 = parse_md(first.optional_field('MD'))
        mism_2 = parse_md(second.optional_field('MD'))
        coord_mism = set()
        for mism in mism_1:
            offset = mism[0]
            coord = first.iv.start + offset
            coord_mism.add(coord)    
        for mism in mism_2:
            offset = mism[0]
            coord = second.iv.start + offset
            coord_mism.add(coord)
                                
        for coord in coord_mism:
            pos = HTSeq.GenomicPosition(genome.id,
                                        coord, strand = first.iv.strand)
            mis[pos] += 1

    # Write coverage and mismatch profiles to file using pickle
    cvg_pickle_fn = sam_fn.replace('.sam', '_cvg.pickle')
    write_pickle(cvg, cvg_pickle_fn)

    mis_pickle_fn = sam_fn.replace('.sam', '_mis.pickle')
    write_pickle(mis, mis_pickle_fn)

    # Create Artemis profile of coverage
    cvg_artemis_fn = sam_fn.replace('.sam', '_cvg.artemis')
    write_artemis(cvg, cvg_artemis_fn)

    mis_artemis_fn = sam_fn.replace('.sam', '_mis.artemis')
    write_artemis(mis, mis_artemis_fn)

    # Create .bedgraph profiles of coverage at plus and minus strands
    cvg_bed_plus_fn = sam_fn.replace('.sam', '_cvg_plus.bedgraph')
    cvg_bed_minus_fn = sam_fn.replace('.sam', '_cvg_minus.bedgraph')
    write_bed(cvg, cvg_bed_plus_fn, cvg_bed_minus_fn)

    mis_bed_plus_fn = sam_fn.replace('.sam', '_mis_plus.bedgraph')
    mis_bed_minus_fn = sam_fn.replace('.sam', '_mis_minus.bedgraph')
    write_bed(mis, mis_bed_plus_fn, mis_bed_minus_fn)

class Feature:
    ''' Contains information about a genomic feature.

Use 'getData' method to collect DMS signal and other information about the feature.
'''
    def __init__(self, sample_name, genome, cvg, mis,
                 feature_type, start, end, strand, locus, name, extend_utr):
        ''' Constructs Feature object.
'''
        
        self.sample_name = sample_name
        self.genome_id = genome.id
        self.type = feature_type
        self.start = int(start) - 1 
        self.end = int(end)
        self.strand = strand
        if self.type == 'five_prime_UTR': #add nucleotides to the 5'UTR to include the first codons
            if self.strand == '+':
                self.end += extend_utr
            else:
                self.start -= extend_utr
        self.locus = locus
        self.name = name
        self.iv = HTSeq.GenomicInterval(self.genome_id, self.start, self.end, self.strand)
        self.coord = range(self.start, self.end)
        self.seq = genome.seq[self.start : self.end]
        self.cvg = list(cvg[self.iv])
        self.mis = list(mis[self.iv])
        self.dms = [None] * len(self.cvg)
        
        if strand == '-':
            self.coord.reverse()
            self.seq = self.seq.reverse_complement()
            self.cvg.reverse()
            self.mis.reverse()

    def __repr__(self):
        return '<Feature object "{0}" from sample "{1}": {2}>'.format(self.locus,
                                            self.sample_name, self.iv)
    
    def getData(self, filter_GT=False):
        ''' Returns a list of lists with the following information:

1) 0-based coordinate of the position.
2) Nucleotide at that coordinate.
3) Coverage at that coordinate.
4) Number of mismatches.
5) DMS signal, if calculated.
'''



        self.data = [list(a) for a in zip(self.coord, list(self.seq), self.cvg, self.mis, self.dms)]
        if filter_GT == True:
            for pos in self.data:
                if pos[1] == 'G' or pos[1] == 'T':
                    pos[3] = 0
        return self.data
       
                    
class Sample:
    '''Contains information about sample.

After DMS profile is calculated, it can be written to .bedtool files with 'writeDMSprofile'
method.
'''
    def __init__(self, sample_name, fasta_fn='NC_003210.1.fa', gff_fn='NC_003210.1.gff3',
                 min_coverage=250, extend_utr=30):
        ''' Generates  Sample object contatining Feature objects with coverage and mismatch data.

Requires:
1) Coverage and mismatch profiles as GenomicArray objects in .pickle files.
The filenames should be: "XXX_cvg.pickle" and "XXX_mis.pickle", where
XXX - sample_name. To generate such files, use 'Generate_profiles.py' script.
2) Genome sequence as .fasta
3) Genomic features in .gff3 format. The attributes field should start
with "ID=XXXX;Name=YYYY;" section.
'''
        self.sample_name = sample_name
        self.dms_profile = None
        self.min_coverage = min_coverage
        self.extend_utr = extend_utr
        # read the genome sequence from fasta file
        self.genome = SeqIO.read(fasta_fn, "fasta")
        self.genome_length = len(self.genome.seq)
        #load GenomicArray objects for coverage and mismatches from .pickle files
        self._cvg_fn = sample_name + '_cvg.pickle'
        self._mis_fn = sample_name + '_mis.pickle'
        with open(self._cvg_fn) as self._CVG:
            self.cvg = pickle.load(self._CVG)
        with open(self._mis_fn) as self._MIS:
            self.mis = pickle.load(self._MIS)
        #read genomic features from gff3 file, create Feature objects, and add these objects to dictionary
        self.features = dict()
        with open(gff_fn) as self._GFF:
            for line in self._GFF:
                line = line.rstrip()
                genome_id, source, feature_type, start, end, score, strand, phase, attributes = line.rsplit('\t')
                locus, name = attributes.rsplit(';')
                locus = locus.replace('ID=', '')
                name = name.replace('Name=', '')
                ft = Feature(self.sample_name, self.genome, self.cvg, self.mis,
                            feature_type, start, end, strand, locus, name, extend_utr=self.extend_utr)
                self.features[ft.locus] = ft
        #determine all feature_types present in sample
        self.feature_types = set()
        for key in self.features.keys():
            ft = self.features[key]
            self.feature_types.add(ft.type)
        #calculate coefficient of DMS reactivity
        self._react = []
        for key in self.features.keys():
            ft = self.features[key]
            if ft.type in ['CDS', 'rRNA', 'tRNA']:
                continue
            data = ft.getData(filter_GT=False)
            for pos in data:
                if pos[2] < self.min_coverage:
                    continue
                if pos[1] == 'A' or pos[1] == 'C':
                    self._react.append(float(pos[3]) / pos[2])
        self.average_dms_signal = numpy.mean(self._react)
        #calculate DMS reactivity values
        for key in self.features.keys():
            ft = self.features[key]
            data = ft.getData()
            ft.dms =[]
            for pos in data:
                if pos[2] < self.min_coverage:
                    ft.dms.append(0)
                elif pos[1] == 'G' or pos[1] == 'T':
                    ft.dms.append(0)
                else:
                    ft.dms.append( (float(pos[3]) / pos[2]) / self.average_dms_signal) 

    def __repr__(self):
        return '<Sample object "{0}">'.format(self.sample_name)

    def calcNucleotideReactivity(self):
        #calculate DMS reactivity for each nucleotide
        nt_react = dict()
        react_A = []
        react_C = []
        react_G = []
        react_T = []
        react = []

        for key in self.features.keys():
            ft = self.features[key]
            data = ft.getData(filter_GT=False)
            for pos in data:
                if pos[2] < self.min_coverage:
                    continue
                if pos[1] == 'A':
                    react_A.append(float(pos[3]) / pos[2])
                elif pos[1] == 'C':
                    react_C.append(float(pos[3]) / pos[2])
                elif pos[1] == 'G':
                    react_G.append(float(pos[3]) / pos[2])
                elif pos[1] == 'T':
                    react_T.append(float(pos[3]) / pos[2])
                if pos[1] == 'A' or pos[1] == 'C':
                    react.append(float(pos[3]) / pos[2])
                    
        nt_react = {'Sample': self.sample_name,
                    'A': numpy.mean(react_A),
                    'C': numpy.mean(react_C),
                    'G': numpy.mean(react_G),
                    'T': numpy.mean(react_T),
                    'AC': numpy.mean(react)}
        
        print 'Coefficients of DMS reactivity in sample', self.sample_name
        print 'Reactivity of A:', nt_react['A']
        print 'Reactivity of C:', nt_react['C']
        print 'Reactivity of G:', nt_react['G']
        print 'Reactivity of T:', nt_react['T']
        print 'Total reactivity (A and C):', nt_react['AC']
        print

        return nt_react

    def cvgByType(self):
        ''' Calculates relative coverage of feature types.
'''
        type_coord = dict()
        for ft_type in self.feature_types:
            for key in self.features.keys():
                ft = self.features[key]
                if ft.type == ft_type:
                    data = ft.getData()
                    for pos in data:
                        type_coord[(pos[0], ft.strand)] = (ft.type, pos[2])
        type_cvg = {}
        for ft_type in self.feature_types:
            type_cvg[ft_type] = 0
        for key in type_coord.keys():
            ft_type, coverage = type_coord[key]
            type_cvg[ft_type] += coverage

        total_cvg = sum(self.cvg.chrom_vectors[self.genome.id]['+']) + sum(self.cvg.chrom_vectors[self.genome.id]['-'])
        print 'Coverage of different feature types in sample {}:'.format(self.sample_name)
        
        for key in type_cvg.keys():
            type_cvg[key] = type_cvg[key] / float(total_cvg)
            print key, '-> {0:.5f}'.format(type_cvg[key])
        print
        
        return type_cvg

    def getFeature(self, name):
        if name in self.features.keys():
            return self.features[name]
        else:
            print 'No feature with the name:', name
    
    def delFeatureType(self, feature_type):
        for key in self.features.keys():
            if self.features[key].type == feature_type:
                del self.features[key]

    def makeFASTA(self):
        ''' Creates .fasta file with sequences of all features except CDS.
'''
        fn =  '{0}_features.fa'.format(self.sample_name)
        FN = open(fn, 'w')
        for key in sorted(self.features.keys()):
            ft = self.features[key]
            FN.write('>{0}_{1}\n'.format(ft.name, self.sample_name))
            FN.write(str(ft.seq)+ '\n')
        FN.close()
    
    def makeSHAPE(self):
        ''' Creates .shape file for all features.
'''
        fn =  '{0}_features.shape'.format(self.sample_name)
        FN = open(fn, 'w')
        for key in sorted(self.features.keys()):
            ft = self.features[key]
            FN.write('>{0}_{1}\n'.format(ft.name, self.sample_name))
            pos = range(1, len(ft.dms)+1)
            data = zip(pos, list(ft.seq), ft.dms, ft.cvg)
            for line in data:
                line = list(line)
                line_cvg = line.pop()
                if line_cvg < self.min_coverage:
                    continue
                if line[1] == 'A' or line[1] == 'C':
                    FN.write('{0}\t{1}\t{2:.5f}\n'.format(*line))
        FN.close()
        
    def makeCOLORMAP(self):
        ''' Creates .colormap file for all features.
'''
        fn =  '{0}_features.colormap'.format(self.sample_name)
        FN = open(fn, 'w')
        for key in sorted(self.features.keys()):
            ft = self.features[key]
            FN.write('>{0}_{1}\n'.format(ft.name, self.sample_name))
            pos = range(1, len(ft.dms)+1)
            data = zip(pos, list(ft.seq), ft.dms, ft.cvg)
            for line in data:
                line = list(line)
                line_cvg = line.pop()
                if line_cvg < self.min_coverage or line[1] == 'T' or line[1] == 'G':
                    line[2] = 50
                FN.write('{0}\t{1}\t{2:.5f}\n'.format(*line))
        FN.close()

    def writeDMSprofile(self):
        self.dms = HTSeq.GenomicArray({self.genome.id: self.genome_length}, stranded=True,
                             storage='ndarray', typecode="d")
        for key in self.features.keys():
            ft = self.features[key]
            for pos in ft.getData():
                self.dms[HTSeq.GenomicPosition(self.genome.id, pos[0], strand=ft.strand)] = pos[4]
        
        dms_bed_plus_fn = self.sample_name + '_dms_plus.bedgraph'
        dms_bed_minus_fn = self.sample_name + '_dms_minus.bedgraph'
        write_bed(self.dms, dms_bed_plus_fn, dms_bed_minus_fn)


class Experiment:

    def __init__(self, name_list, fasta_fn='NC_003210.1.fa', gff_fn='NC_003210.1.gff3',
                 min_coverage=250, extend_utr=30):
        self.name_list = name_list
        self.fasta_fn = fasta_fn
        self.gff_fn = gff_fn
        self.min_coverage = min_coverage
        self.extend_utr = extend_utr

    def calcStat(self):
        REACT = open('Nucleotide_reactivity.txt', 'w')
        REACT.write('Sample\tA\tC\tG\tT\tAC\n')
        TYPECVG = open('Coverage_by_type.txt', 'w')
        TYPECVG.write('Sample\tSRP_RNA\trRNA\ttRNA\tsmall_regulatory_ncRNA\tfive_prime_UTR\tantisense_RNA\tCDS\tRNase_P_RNA\tribosome_entry_site\tRiboswitches_and_attenuators\n')
        for name in self.name_list:
            print name
            s = Sample(name, self.fasta_fn, self.gff_fn)
            type_cvg = s.cvgByType()
            type_cvg['Sample'] = s.sample_name
            TYPECVG.write('{Sample}\t{SRP_RNA}\t{rRNA}\t{tRNA}\t{small_regulatory_ncRNA}\t{five_prime_UTR}\t{antisense_RNA}\t{CDS}\t{RNase_P_RNA}\t{ribosome_entry_site}\t{transcription_regulatory_region}\n'.format(**type_cvg))
            nt_react = s.calcNucleotideReactivity()
            REACT.write('{Sample}\t{A}\t{C}\t{G}\t{T}\t{AC}\n'.format(**nt_react))
        REACT.close()
        TYPECVG.close()
        
    def exportData(self):
        '''Exports data to files.
1) .shape file suitable for RNAfold program
2) .colormap file suitable for VARNA program
3) .bedgraph files with DMS reactivity profiles
4) 'Experiment.txt' file with coverage and mismatch data for further statistical analysis.
'''
        self.header = ['Position_genome', 'Position_TSS', 'Position_start_codon',
                       'Strand', 'Feature', 'Type', 'Nucleotide']
        self.exp_dict = dict()
        data_fn = 'Experiment.txt'
        
        for name in self.name_list:
            print name
            s = Sample(name, self.fasta_fn, self.gff_fn,
                       min_coverage=self.min_coverage, extend_utr=self.extend_utr)
            
            s.makeSHAPE()
            s.makeCOLORMAP()
            s.writeDMSprofile()
            
            self.header.append('cvg_' + s.sample_name.replace('lib_', ''))
            self.header.append('mis_' + s.sample_name.replace('lib_', ''))
            sample_dict = dict()
            
            if name == self.name_list[0]:
                for key in s.features.keys():
                    ft = s.features[key]
                    if ft.type in ['CDS', 'rRNA', 'tRNA', 'ribosome_entry_site']:
                        continue
                    
                    for pos in ft.getData():
                        coord, nt, cvg, mis, dms = pos
                        coord += 1
                        
                        if ft.strand == '+':
                            tss_rel = str(coord - ft.start)
                            start_rel = str(coord - ft.end + self.extend_utr)
                        else:
                            tss_rel = str(ft.end - coord + 1)
                            start_rel = str((ft.start + 1) - coord + self.extend_utr)
                        if ft.type != 'five_prime_UTR':
                            start_rel = 'NA'

                        if nt == 'G' or nt == 'T':
                            continue
                        if coord in sample_dict.keys():
                            sample_dict[coord][4] = sample_dict[coord][4] + ', {0}({1})'.format(ft.locus, ft.name)
                            sample_dict[coord][5] = sample_dict[coord][5] + ', ' + ft.type
                            if ft.type == 'five_prime_UTR':
                                sample_dict[coord][1] = tss_rel
                                sample_dict[coord][2] = start_rel
                        else:
                            sample_dict[coord] = [str(coord), tss_rel, start_rel, ft.strand, 
                                                  '{0}({1})'.format(ft.locus, ft.name),
                                                  ft.type, nt, str(cvg), str(mis)]
                            
                for sample_key in sample_dict.keys():
                    self.exp_dict[sample_key] = sample_dict[sample_key]
            else:
                for key in s.features.keys():
                    ft = s.features[key]
                    if ft.type in ['CDS', 'rRNA', 'tRNA', 'ribosome_entry_site']:
                        continue
                    for pos in ft.getData():
                        coord, nt, cvg, mis, dms = pos
                        coord += 1
                        if nt == 'G' or nt == 'T':
                            continue
                        if coord in sample_dict.keys():
                            continue
                        else:
                            sample_dict[coord] = [str(cvg), str(mis)]
                for sample_key in sample_dict.keys():
                    self.exp_dict[sample_key].extend(sample_dict[sample_key])
        
        
        print 'Writing data to file...'
        with open(data_fn, 'w') as FN:
            FN.write('\t'.join(self.header) + '\n')
            for key in sorted(self.exp_dict.keys()):
                FN.write('\t'.join(self.exp_dict[key]) + '\n')

Exp = Experiment(['lib_K', 'lib_invitro', 'lib_total',
		  'lib_37A', 'lib_37B',
		  'lib_26', 'lib_26shift', 'lib_37shift',  
		    'lib_hfq',  'lib_lhrA',  'lib_prfA'])

Exp_60 = Experiment(['lib_K', 'lib_invitro', 'lib_37A', 'lib_37B'], extend_utr=60)
