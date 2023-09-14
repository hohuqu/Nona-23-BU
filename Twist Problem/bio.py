from collections import defaultdict
import re
import itertools
from typing import NamedTuple
import pandas as pd
from Bio.Seq import Seq
import Bio.Restriction as rst


class Fragment(NamedTuple):
    seq: Seq
    end5: rst.RestrictionType
    end3: rst.RestrictionType

    def can_ligate(self, other):
        compat = self.end3.compat % other.end5.compat
        overhang_zero = self.overhang + other.overhang == 0
        return compat and overhang_zero

    def ligate(self, other):
        if not self.can_ligate(other):
            raise ValueError('Cannot ligate')
        return Fragment(self.seq + other.seq, self.end5, other.end3)


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


def equivalences(stuff, equiv_rel):
    equivalences = []
    for e in stuff:
        found = False
        for equiv in equivalences:
            if equiv_rel(e, equiv[0]):
                equiv.append(e)
                found = True
        if not found:
            equivalences.append([e])
    return {e[0]: e for i, e in enumerate(equivalences)}


def chew_seq(seq, enzymes, frags):
    prev_e = None
    blunt_example = rst.BsrBI
    for i, e in enumerate(enzymes):
        s1, s2 = e.catalyze(seq)
        seq = s2
        if i == 0:
            e1 = rst.BsrBI
            e2 = e
        elif i == len(enzymes) - 1:
            e1 = prev_e
            e2 = rst.BsrBI)
        else:
            e1 = prev_e
            e2 = End(..., e)
        prev_e = e
        frags.append(Fragment(Seq(s1), e1, e2))
    return frags


def digestion_ligation(seqs, enzymes_2p, equiv_2p):
    analyses = {i: rst.Analysis(enzymes_2p, s) for i, s in enumerate(seqs)}
    equiv_classes = defaultdict(set)
    for i, a in analyses.items():
        for e in a.with_sites():
            k = next((k for k in equiv_2p if e % k), None)
            if k:
                equiv_classes[i].add(k)
    breakpoint()
    assert all(len(s) == 1 for s in equiv_classes.values()), equiv_classes


def main():
    # seq = Seq('AAAAGAATTCAAAAAAAA')
    # mat = enzyme_equiv_matrix()
    enzymes_2p = [e for e in rst.AllEnzymes if not is2s(e)]
    enzymes_2s = [e for e in rst.AllEnzymes if is2s(e)]
    isoclass = equivalences(enzymes_2p, lambda e1, e2: (e1 % e2) and e1.is_equischizomer(e2))
    compatclass = equivalences(isoclass.keys(), lambda e1, e2: e1 % e2)
    breakpoint()
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
