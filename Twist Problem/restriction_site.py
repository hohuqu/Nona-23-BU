from collections import defaultdict
import re
import itertools
from typing import NamedTuple

from Bio.Seq import Seq
import Bio.Restriction as rst
from tqdm import tqdm


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


BLUNT_EXAMPLE = rst.BsrBI
ENZYMES_2P = [e for e in rst.AllEnzymes if not is2s(e)]
ENZYMES_2S = [e for e in rst.AllEnzymes if is2s(e)]
ISOCLASS = equivalences(ENZYMES_2P, lambda e1, e2: (e1 % e2) and e1.is_equischizomer(e2))


class Fragment(NamedTuple):
    seq: Seq
    end5: object  # enzyme
    end3: object  # enzyme

    def can_ligate(self, other):
        compat = self.end3.compat % other.end5.compat
        overhang_zero = self.overhang + other.overhang == 0
        return compat and overhang_zero

    def ligate(self, other):
        if not self.can_ligate(other):
            raise ValueError('Cannot ligate')
        return Fragment(self.seq + other.seq, self.end5, other.end3)


def chew_seq(seq, enzymes, frags):
    BLUNT_EXAMPLE = rst.BsrBI
    prev_e = BLUNT_EXAMPLE
    for i, e in enumerate(enzymes):
        try:
            (s1, *s2) = e.catalyze(seq)
        except NotImplementedError:
            continue
        s2 = Seq(''.join(map(str, s2)))
        if i == 0:
            e1 = BLUNT_EXAMPLE
            e2 = e
        elif i == len(enzymes) - 1:
            e1 = e
            e2 = BLUNT_EXAMPLE
        else:
            e1 = prev_e
            e2 = e
        prev_e = e
        frags.add(Fragment(Seq(s1), e1, e2))
    return frags


def with_sites_inverse(d):
    res = defaultdict(list)
    for e, sites in d.items():
        for s in sites:
            res[s].append(e)
    return res


def sites_to_equiv_reps(a, equiv_2p):
    res = {}
    for e, sites in a.with_sites().items():
        k = next((k for k in equiv_2p if e % k), None)
        if k:
            res[k] = sites
    return res


def digestion_ligation(seqs, enzymes_2p, equiv_2p):
    frags = set()
    for s in tqdm(seqs):
        a = rst.Analysis(enzymes_2p, s)
        sites = sites_to_equiv_reps(a, equiv_2p)
        inverse = with_sites_inverse(sites)
        for es in tqdm(itertools.product(*inverse.values())):
            chew_seq(s, es, frags)
    return frags


def build(frags, max_recur):
    seqs = []

    def _build(seq, end5, recur=1):
        if end5 is BLUNT_EXAMPLE or recur == max_recur:
            seqs.append(seq)
            return
        for f in (f for f in frags if f.end5 % end5):
            _build(seq + f.seq, f.end3, recur=recur + 1)

    for f in (f for f in frags if f.end5 is BLUNT_EXAMPLE):
        _build(f.seq, f.end3)
    return seqs


def restriction(sequence1, sequence2, max_recur=5):
    frags = digestion_ligation([sequence1, sequence2], ENZYMES_2P, ISOCLASS)
    return build(frags, max_recur)


def main():
    seq1 = Seq('TCCCTGGGCTCTTTTAGTGGACGGAGACCCAGCTGTCAGTTTGTTGTAATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA')
    seq2 = Seq('CTGCCCAAGCCTACCGTGAATCATCTAATCCCTCCATGGAGTAAGTGGTGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT')
    # seq1 = Seq('AAAAAGATCCAAAAAA')
    # seq2 = Seq('TTTTGGATCCTTTTTT')
    print(len(restriction(seq1, seq2)))


if __name__ == '__main__':
    main()
