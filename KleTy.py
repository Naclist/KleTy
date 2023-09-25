#! /usr/bin/python3

import click, tempfile, gzip, os, sys


from amrsearch import amrsearch
from cgMLST import cgmlst
from PlasT import plast

def readFile(fname) :
    fin = gzip.open(fname, 'rt') if fname.upper().endswith('.GZ') else open(fname, 'rt')
    text = fin.readlines()
    fin.close()
    if text[0].startswith('@') :
        text = [ (t if i % 4 == 1 else '>' + t[1:]) for i, t in enumerate(text) if i % 4 < 2 ]
    return ''.join(text)


def write_cgMLST(prefix, alleles) :
    fout = '{0}.cgMLST.profile.gz'.format(prefix)
    with gzip.open(fout, 'wt') as fout :
        fout.write('#Query\t{0}\n'.format('\t'.join([g for g in sorted(alleles.keys())])))
        fout.write('{0}\t{1}\n'.format(prefix, 
                                       '\t'.join([allele.get('value_md5', '-') for g, allele in sorted(alleles.items())])))





@click.command()
@click.option('-q', '--query', help='query genome in fasta or fastq format. May be gzipped.', required=True)
@click.option('-o', '--prefix', help='prefix for output. default: query filename', default=None)
@click.option('-n', '--n_proc', help='number of process to use. default: 8', default=8, type=int)
@click.option('-m', '--skip_gene', help='flag to skip AMR/VF searching. default: False', default=False, is_flag=True)
@click.option('-g', '--skip_cgmlst', help='flag to skip cgMLST. default: False', default=False, is_flag=True)
@click.option('-p', '--skip_plasmid', help='flag to skip plasmid typing. default: False', default=False, is_flag=True)
def klebtyper(query, prefix, skip_gene, skip_cgmlst, skip_plasmid, n_proc) :
    if prefix == None :
        prefix = os.path.basename(query).rsplit('.', 1)[0]
    
    with tempfile.TemporaryDirectory(dir='.', prefix='kt_') as tmpdir :
        tmpfile = os.path.join(tmpdir, 'qry.fna')
        with open(tmpfile, 'wt') as fout :
            fout.write(readFile(query))

        if not skip_plasmid :
            sys.stderr.write('Running plasmid...')
            plasmids = plast(tmpfile, n_proc=n_proc)
            plasmids = [p for p in plasmids if p[0].startswith('PT')]
            sys.stderr.write('Done.\n')
        else :
            plasmids = []
        
        if not skip_gene :
            sys.stderr.write('Running AMR/VF gene searching...')
            genes = amrsearch(tmpfile, n_proc=n_proc)
            sys.stderr.write('Done.\n')
        else :
            genes = []
        
        cont_replicon = {cont:i for i, plasmid in enumerate(plasmids) for cont in plasmid[3].split(',') }

        replicons = {ri:{} for ri in range(-1, len(plasmids)) }
        for gene in genes :
            if gene[8] in ('oqxA', 'oqxB', 'oqxR', 'fosA5_fam', 'ramR', 'uhpT') or gene[9] in ('pmrB_R256G'):
                continue
            #elif gene[8] in ('blaSHV', 'blaLEN', 'blaOKP') and 'ESBL' not in gene[9] :
             #   continue
            
            replicon = cont_replicon.get(gene[0], len(plasmids))
            for category in set(gene[9].split('/') + gene[10].split('/')) :
                if category not in {'', 'MULTIDRUG', 'ORGANOMERCURY', 'PHENYLMERCURY', 'POLYKETIDE', 'ARSENITE', 'ARSENATE'} :
                    if category in {'ARSENIC', 'COPPER', 'SILVER', 'GOLD', 'FLUORIDE', 'MERCURY', 'TELLURIUM', 'NICKEL', 'QUATERNARY_AMMONIUM'} :
                        category = 'STRESS:{0}'.format(category)
                    elif category in ('INC_TYPE', 'MOB_TYPE', 'MPF_TYPE') :
                        category = 'REPLICON:{0}'.format(category)
                    elif category in ('iuc', 'rmp', 'ybt', 'clb', 'iro') :
                        category = 'VIRULENCE:{0}'.format(category)
                    else :
                        category = 'AMR:{0}'.format(category)
                    
                    if replicon not in replicons :
                        replicons[replicon] = {}
                    if category not in replicons[replicon] :
                        replicons[replicon][category] = {gene[7]}
                    else :
                        replicons[replicon][category].update({gene[7]})
                    if category not in replicons[-1] :
                        replicons[-1][category] = {gene[7]}
                    else :
                        replicons[-1][category].update({gene[7]})

        if skip_plasmid :
            replicons.pop(0, None)

        resistances = 'AMR:AMINOGLYCOSIDE|AMR:AMIKACIN|AMR:APRAMYCIN|AMR:GENTAMICIN|AMR:HYGROMYCIN|AMR:KANAMYCIN|AMR:SPECTINOMYCIN|AMR:STREPTOMYCIN|AMR:TOBRAMYCIN|AMR:BETA-LACTAM|AMR:CARBAPENEM|AMR:CEPHALOSPORIN|AMR:ESBL|AMR:INHIBITOR-RESISTANT|AMR:COLISTIN|AMR:FOSFOMYCIN|AMR:LINCOSAMIDE|AMR:MACROLIDE|AMR:PHENICOL|AMR:CHLORAMPHENICOL|AMR:FLORFENICOL|AMR:KIRROMYCIN|AMR:PULVOMYCIN|AMR:QUINOLONE|AMR:RIFAMYCIN|AMR:STREPTOTHRICIN|AMR:SULFONAMIDE|AMR:TETRACYCLINE|AMR:TIGECYCLINE|AMR:TRIMETHOPRIM|AMR:BLEOMYCIN|STRESS:COPPER|STRESS:MERCURY|STRESS:NICKEL|STRESS:SILVER|STRESS:TELLURIUM|STRESS:ARSENIC|STRESS:FLUORIDE|STRESS:QUATERNARY_AMMONIUM|VIRULENCE:clb|VIRULENCE:iro|VIRULENCE:iuc|VIRULENCE:rmp|VIRULENCE:ybt|Others|REPLICON:INC_TYPE|REPLICON:MOB_TYPE|REPLICON:MPF_TYPE'.split('|')
        for r, drugs in replicons.items() :
            for drug in resistances :
                replicons[r][drug] = ','.join(sorted(replicons[r][drug])) if drug in replicons[r] else '-'
            for drug in sorted(drugs.keys()) :
                if drug not in resistances :
                    x = '{0}({1})'.format(','.join(sorted(replicons[r][drug])), drug)
                    if replicons[r]['Others'] == '-' :
                        replicons[r]['Others'] = x
                    else :
                        replicons[r]['Others'] += '|' + x
        
        if not skip_cgmlst :
            sys.stderr.write('Running cgMLST...')
            alleles, hiercc = cgmlst(tmpfile, n_thread=n_proc)
            write_cgMLST(prefix, alleles)
            sys.stderr.write('Done.\n')
        else :
            hiercc = {}

        with open(prefix + '.KleTy', 'wt') as fout :
            fields = ['REPLICON', 'SPECIES', 'HC1360.500.200.100.50.20.10.5.2', 'REFERENCE', 'PLASTYPE', 'COVERAGE'] + \
                resistances + ['ANNOTATION', 'CONTIGS']
            fout.write('\t'.join(fields)+'\n')
            for rep_id in sorted(replicons.keys()) :
                if rep_id == -1 :
                    res = [prefix+':ALL', hiercc.get('species', 'ND'), hiercc.get('HC1360.500.200.100.50.20.10.5.2', 'ND'), hiercc.get('reference','ND'), '-', '-'] + \
                        [replicons[-1][r] for r in resistances] + ['-', '-']
                elif rep_id >= len(plasmids) :
                    res = [prefix+':Others', '-', '-', '-', '-', '-'] + [replicons[rep_id][r] for r in resistances] + \
                        ['-', '-']
                else :
                    plasmid = plasmids[rep_id]
                    ref, ann = plasmid[-1].split(' ', 1)
                    ref = ref.split('_', 2)[2]
                    res = ['{0}:P{1}'.format(prefix, rep_id+1), '-', '-', ref, plasmid[0].replace('\t', ','), plasmid[1]] + \
                        [replicons[rep_id][r] for r in resistances] + [ann.replace(' ', '_'), plasmid[3]]
                fout.write('\t'.join(res)+'\n')


if __name__ == '__main__' :
    klebtyper()
