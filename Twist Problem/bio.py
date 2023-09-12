from collections import defaultdict
import re
import itertools
import pandas as pd
from Bio.Seq import Seq
import Bio.Restriction as rst


def is2s(e):
    """
    Is Type 2 enzyme
    """
    s = e.elucidate()
    find = re.search(r'[\^_N]', s)
    if find:
        i = find.span(0)[0]
        return any(set(sp).issubset({'^', '_', 'N'}) for sp in (s[:i], s[i:]))
    return False


def enzyme_equivalences(enzymes):
    equivalences = []
    for e in enzymes:
        found = False
        for equiv in equivalences:
            if e % equiv[0]:
                equiv.append(e)
                found = True
        if not found:
            equivalences.append([e])
    return {str(e[0]): e for i, e in enumerate(equivalences)}


def enzyme_equiv_matrix():
    # equivs = [e1 % e2 for e1, e2 in itertools.product(enzymes, enzymes)]
    enzyme_names = rst.AllEnzymes.elements()
    df = pd.DataFrame(columns=enzyme_names, index=enzyme_names, dtype=bool)
    df[:] = False
    for e1, e2 in itertools.product(enzyme_names, enzyme_names):
        if rst.AllEnzymes.get(e1) % rst.AllEnzymes.get(e2):
            df[e1].loc[e2] = True
    return df


def digestion_ligation(seqs, enzymes_2p, equiv_2p):
    analyses = {i: rst.Analysis(enzymes_2p, s) for i, s in enumerate(seqs)}
    equiv_classes = defaultdict(list)
    for i, a in analyses.items():
        for e in a.with_sites().keys():
            for k in (rst.AllEnzymes.get(k) for k in equiv_2p):
                if e % k:
                    equiv_classes[i].append(str(e))
    breakpoint()
    assert all(len(s) == 1 for s in equiv_classes.values()), equiv_classes


def main():
    # seq = Seq('AAAAGAATTCAAAAAAAA')
    # mat = enzyme_equiv_matrix()
    enzymes_2p = [e for e in rst.AllEnzymes if not is2s(e)]
    enzymes_2s = [e for e in rst.AllEnzymes if is2s(e)]
    equiv_2p = enzyme_equivalences(enzymes_2p)
    equiv_2s = enzyme_equivalences(enzymes_2s)
    seq1 = Seq('AAAAAGATCCAAAAAA')
    seq2 = Seq('TTTTGGATCCTTTTTT')
    digestion_ligation([seq1, seq2], enzymes_2p, equiv_2p)
    print(list(rst.AllEnzymes))
    # print(f'{rst.XhoII.elucidate()=}')
    # print(f'{rst.XhoII.elucidate()=}')
    # res = rst.Analysis(rst.AllEnzymes, seq)
    # res.print_that()
    print(f'{(rst.BsaI % rst.XhoII)=}')
    # print(f'{(rst.BsaI % rst.BcoDI)=}')


if __name__ == '__main__':
    main()
