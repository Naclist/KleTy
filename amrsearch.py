import os, click, pandas as pd, numpy as np, re
from Bio import Seq
try :
    from uberBlast import uberBlast, readFastq
    from configure import dbs
except :
    from .uberBlast import uberBlast, readFastq
    from .configure import dbs


resfams = dbs['resfams']
dbs = dbs['amr_dbs']


def get_tops(bsn) :
    bsn = bsn[np.argsort(-bsn.T[11]*bsn.T[3]/bsn.T[12])]
    covs = { c:[] for c in np.unique(bsn.T[1]) }
    for p in bsn :
        s, e = p[8:10] if p[8] < p[9] else (p[9], p[8])
        x = np.zeros(e-s+1, dtype=np.uint8)
        for c in covs[p[1]] :
            if c[0] > e or c[1] < s :
                continue
            ovl = [max(c[0], s), min(c[1], e)]
            if ovl[1] - ovl[0] + 1 >= 0.3 * (c[1] - c[0] + 1) :
                p[1] = ''
                break
            x[ovl[0]-s:ovl[1]-s+1] = 1
            if np.sum(x) >= 0.3 * (e - s + 1) :
                p[1] = ''
                break
        else :
            covs[p[1]].append([s, e])
    return bsn[bsn.T[1] != '']
    

def rc(s) :
    c = {'A':'T', 'T':'A', 'G':'C', 'C':'G', '-':'-'}
    return ''.join([ c.get(x, 'N') for x in s[::-1].upper() ])


def amrsearch(query, g_table=11, n_proc=8) :
    bla_category = {}
    if os.path.isfile(resfams) :
        with open(resfams, 'rt') as fin :
            for line in fin :
                p = line.strip().split('|')
                variant = p[4]
                if variant in bla_category :
                    continue
                func, t = [], p[7]
                if 'carbapenem' in t :
                    func.append('CARBAPENEM')
                if 'extended-spectrum' in t :
                    func.append('ESBL')
                if 'inhibitor-resistant' in t or 'metallo-' in t :
                    func.append('INHIBITOR-RESISTANT')
                if len(func) :
                    bla_category[variant] = func
    
    qry, _ = readFastq(query)
    res = []
    resistance = {}

    for title, (db, site, min_iden, min_cov) in dbs.items() :
        if title.startswith('prot') :
            bsn = uberBlast('-r {0} -q {1} -f --diamondx --min_id {2} --min_ratio {3} -t {4} -p -s 2 -e 0,3 -m --merge_gap 300'.format(
                query, db, min_iden/100., min_cov/100., n_proc).split())
            if title == 'prot' :
                for p in bsn :
                    ref = p[0].split(' ')[0].split('|')
                    if ref[6] != 'mutation' :
                        resistance[ref[4]] = ref[8:10]
                        resistance[ref[5]] = ref[8:10]

        else :
            bsn = uberBlast('-r {0} -q {1} -f --blastn --min_id {2} --min_ratio {3} -t {4} -p -s 2 -e 0,3 -m --merge_gap 300'.format(
                query, db, min_iden/100., min_cov/100., n_proc).split())
            
        bsn = get_tops(bsn)
        
        if site == None :
            for p in bsn:
                if title == 'nucl' :
                    ref = p[0].split(' ')[0].split('|')
                    category = ''
                elif title == 'nuc_vir' :
                    if p[0].find('delete') >= 0 :
                        continue
                    if p[0].find('__rmp') >= 0 :
                        p[0] = re.split('__', p[0])[2]
                    x = p[0].rsplit('_', 1)[0]
                    ref = ['', p[0], '', '', '', p[0], x, x, x]
                    category = p[0][:3]
                    if category in ('fyu', 'irp') :
                        category = 'ybt'
                    elif category in ('iut') :
                        category = 'iuc'
                elif title.endswith('plasmids') :
                    x = p[0].rsplit('|', 1)[-1]
                    ref = ['', p[0], '', '', '', x, x, x, x]
                    category = 'MOB_TYPE' if x.startswith('MOB') else ('MPF_TYPE' if x.startswith('MPF') else 'INC_TYPE')
                else :
                    x = p[0].rsplit('_', 1)[0]
                    ref = ['', p[0], '', '', '', p[0], x, x, x]
                    category = p[0][:3]

                qry_na = Seq.Seq(qry[p[1]][p[8]-1:p[9]]) if p[8] < p[9] else Seq.Seq(qry[p[1]][p[9]-1:p[8]]).reverse_complement()
                if title == 'nuc_plasmids' :
                    res.append([p[1], p[8], p[9], p[2]*100, np.round(100.*(p[7]-p[6]+1)/p[12], 2), ref[1], 'BLASTN', ref[5], ref[6], category, category, ref[7]])
                elif title == 'prot_plasmids' :
                    res.append([p[1], p[8], p[9], p[2]*100, np.round(100.*(p[7]-p[6]+1)/p[12], 2), ref[1], 'BLASTN', ref[5], ref[6], category, category, ref[7]])
                elif (abs(p[9] - p[8]) + 1) % 3 > 0 :
                    if (max([len(x) for x in qry_na[:int((len(qry_na))/3)*3].translate(g_table).split('*')])*3) >= 0.85 * len(qry_na) or \
                        (max([len(x) for x in qry_na[len(qry_na)%3:].translate(g_table).split('*')])*3) >= 0.85 * len(qry_na) :
                        res.append([p[1], p[8], p[9], p[2]*100, np.round(100.*(p[7]-p[6]+1)/p[12], 2), ref[1], 'BLASTN', ref[5], ref[6], category, category, ref[7]])
                    else :
                        res.append([p[1], p[8], p[9], p[2]*100, np.round(100.*(p[7]-p[6]+1)/p[12], 2), ref[1], 'BLASTN', ref[5] + '(*Frameshift)', ref[6], category, category, ref[7]])
                else :
                    qry_aa = qry_na.translate(g_table)
                    if '*' in qry_aa[:-1] and max([len(x) for x in qry_aa.split('*')])*3 < 0.85 * len(qry_na) :
                        res.append([p[1], p[8], p[9], p[2]*100, np.round(100.*(p[7]-p[6]+1)/p[12], 2), ref[1], 'BLASTN', ref[5] + '(*Premature)', ref[6], category, category, ref[7]])
                    else :
                        res.append([p[1], p[8], p[9], p[2]*100, np.round(100.*(p[7]-p[6]+1)/p[12], 2), ref[1], 'BLASTN', ref[5], ref[6], category, category, ref[7]])
        else :
            sites = pd.read_csv(site, sep='\t').values[:, -6:]
            for p in bsn :
                ref = p[0].split(' ')[0].split('|')
                if len(ref) > 1 :
                    if ref[6] != 'mutation' :
                        res.append([p[1], p[8], p[9], p[2]*100, np.round(100.*(p[7]-p[6]+1)/p[12], 2), ref[1], 'BLASTX', ref[4], ref[5], ref[8], ref[9], ref[10]])
                    ref = ref[1:]
                psite = {}
                for x in sites[sites.T[0] == ref[0]] :
                    loc = int(x[1])
                    if loc >= p[6] and loc <= p[7] :
                        variation = list(re.findall(r'_([^-\d]+)[-\d]+([^-\d]+)', x[2])[0])
                        if variation[1].upper() == 'DEL' :
                            variation[1] = '-' * len(variation[0])
                        elif variation[1].upper() == 'STOP' :
                            variation[1] = '*'
                        if len(variation[0]) > len(variation[1]) :
                            variation[1] += '-' * (len(variation[0]) - len(variation[1]))
                        elif len(variation[1]) > len(variation[0]) :
                            variation[0] += '-' * (len(variation[1]) - len(variation[0]))
                        psite[loc] = [x, variation]
                
                if len(psite) < 1 :
                    continue
                qry_s = Seq.Seq(qry[p[1]][p[8]-1:p[9]]) if p[8] < p[9] else Seq.Seq(qry[p[1]][p[9]-1:p[8]]).reverse_complement()
                cigar = [[int(s), t] for s, t in re.findall('(\d+)([MID])', p[14])]
                if title == 'prot' :
                    mv = 3
                    qry_s = qry_s.translate(g_table)
                    cigar = [ [int(s/3), t] for s, t in cigar ]
                else :
                    mv = 1
                c = 0
                qq = []
                for s, t in cigar :
                    if t != 'I' :
                        qq.append(str(qry_s[c:c+s]))
                        c += s
                    else :
                        qq.append('-' * s)
                qq = ''.join(qq)
                ss = ''.join([ (t * s) if t != 'D' else ('-'*s) for s, t in cigar ])
                
                if p[8] > p[9] :
                    mv *= -1
                qi, ri = p[8] - mv, p[6] - 1
                
                for xi, (q, r) in enumerate(zip(list(qq), list(ss))) :
                    if q != '-' :
                        qi += mv
                    if r != '-' :
                        ri += 1
                    if ri in psite :
                        x, variation = psite.pop(ri)
                        bq = qq[xi:xi+len(variation[1])]
                        if bq == variation[1] : 
                            len_v = len(variation[1].replace('-', ''))
                            len_x = 1 if p[8] > p[9] else -1
                            res.append([p[1], qi, qi+mv*len_v+len_x, p[2]*100, np.round(100.*(p[7]-p[6]+1)/p[12], 2), x[0], 'MUTATION', x[2], x[2].split('_', 1)[0], x[4], x[3], x[5]])
                        if len(psite) < 1 :
                            break
        # import datetime
        # dt_string = datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        # print(title, dt_string)

    for r in res :
        if r[9] == '' :
            if r[7] in resistance :
                r[9:11] = resistance.get(r[7], ['', ''])
            else :
                r[9:11] = resistance.get(r[8], ['', ''])
        if r[7] in bla_category :
            r[9] = '/'.join( sorted(set(r[9].split('/')+bla_category[r[7]])) )
    res = sorted([ r for r in res if r[9] != '' ])
    for j, r2 in list(enumerate(res))[1:] :
        todel = []
        for r1 in res[j-1::-1] :
            if r1[0] == '' :
                continue
            if r1[0] != r2[0] :
                break
            s1, e1 = sorted(r1[1:3])
            s2, e2 = sorted(r2[1:3])
            ovl = min(e2, e1) - max(s1, s2) + 1
            if ovl >= 0.5 * (e1-s1+1) or ovl >= 0.5 * (e2-s2+1) :
                sc1 = (e1-s1+1) * r1[3] * r1[4]/10000.
                sc2 = (e2-s2+1) * r2[3] * r2[4]/10000.
                if sc1 >= sc2 :
                    r2[0] = ''
                    break
                else :
                    todel.append(r1) #r1[0] = ''
        if r2[0] != '' :
            for r1 in todel :
                r1[0] = ''
    return [r for r in res if r[0] != '']


@click.command()
@click.option('-q', '--query')
def main(query) :
    amrsearch(query)



if __name__ == '__main__' :
    main()