import sys
from collections import defaultdict, OrderedDict, Counter

__version__ = "0.2.5"


class Seq:
    def __init__(self, tag, seq):
        self.tag = tag
        self.seq = seq


class Tag(OrderedDict):
    def make(self, tag):
        cells = tag.strip().split('|')
        cell_len = len(cells)
        n = 0
        while n + 1 < cell_len:
            self[cells[n]] = cells[n+1]
            n += 2

    def __str__(self):
        tags = []
        for k in self:
            tags.extend([str(k), self[k]])
        return '|'.join(tags)

    def add(self, k, v):
        self[k] = v.replace('|', '_')

    @classmethod
    def parse(cls, tag):
        ob = cls()
        ob.make(tag)
        return ob


class MSeq(list):
    def parse_fasta(self, string):
        tag, seq = '', ''
        for l in string.splitlines():
            if not l:
                continue
            if l[0] == '>':
                if seq:
                    self.add(tag.strip(), seq.strip().replace(' ', ''))
                    seq = ''
                tag = l[1:]
            else:
                seq += l
        self.add(tag.strip(), seq.strip().replace(' ', ''))

    def parse_phylip(self, string):
        blocks = string.split('\n\n')
        i = 0
        for block in blocks:
            if not block:
                continue
            j = 0
            for l in block.splitlines():
                if not l:
                    continue
                if l[0] != ' ':
                    cs = l.split()
                    self.add(cs[0], ''.join(cs[1:]))
                elif i > 0:
                    self[j].seq += l.strip().replace(' ', '')
                j += 1
            i += 1

    def parse_msf(self, string):
        _, alg = string.split('//')
        blocks = alg.strip().split('\n\n')
        tags = []
        tag2seq = {}
        for block in blocks:
            if not block:
                continue
            for l in block.splitlines():
                if not l:
                    continue
                cs = l.split()
                if not cs[0] in tags:
                    tags.append(cs[0])
                    tag2seq[cs[0]] = ''
                tag2seq[cs[0]] += ''.join(cs[1:])
        for tag in tags:
            self.add(tag, tag2seq[tag])


    def iterparse_fasta(self, ITER):
        tag, seq = '', ''
        for l in ITER:
            l = l.strip()
            if l and l[0] == '>':
                if seq:
                    yield Seq(tag.rstrip(), seq.strip().replace(' ', ''))
                    seq = ''
                tag = l[1:]
            else:
                seq += l
        yield Seq(tag.rstrip(), seq.strip().replace(' ', ''))

    def add(self, tag, seq):
        self.append(Seq(tag, seq))
    
    def erase_common_gaps(self):
        pos2num  = defaultdict(int)
        for o in self:
            for idx in [i for i, x in enumerate(o.seq) if x == '-']:
                pos2num[idx] += 1
        seqnum = len(self)
        pos_to_erase = []
        for k, v in pos2num.items():
            if v == seqnum:
                pos_to_erase.append(k)
        pos_to_erase.sort(reverse=True)
        for o in self:
            for pos in pos_to_erase:
                o.seq = o.seq[:pos] + o.seq[pos+1:]

    def upper(self):
        for o in self:
            if not o.seq.isupper():
                o.seq = o.seq.upper()

    def lower(self):
        for o in self:
            if not o.seq.islower():
                o.seq = o.seq.lower()

    def rm(self, char):
        for o in self:
            o.seq = o.seq.replace(char, '') 

    def seq_type(self):
        seq_len = 0
        c = Counter(self[0].seq.lower())
        for k, v in c.items():
            if k not in ('-', 'x'):
                seq_len += v
        if (c['a'] + c['c'] + c['g'] + c['t'] + c['u']) / float(seq_len) > .8:
            if c['u']:
                return 'RNA'
            else:
                return 'DNA'
        else:
            return 'Protein'

    def seq_len(self, raise_error=False):
        seqLen = 0
        for ob in self:
            currentLen = len(ob.seq)
            if seqLen:
                if seqLen != currentLen:
                    if raise_error:
                        sys.stderr.write('Mismatch in sequence len.\n')
                        sys.exit(0)
                    else:
                        return False
            else:
                seqLen = currentLen
        return seqLen

    def write_phylip(
        self, oh,
        max_tag_len=10, line_len=60, block_len=10
    ):
        seq_len = len(self[0].seq)
        seq_num = len(self)
        seq_offset = max_tag_len+3
        oh.write(' %s %s\n' % (seq_num, seq_len))
        for i, s in enumerate(range(0, seq_len, line_len)):
            for j in range(seq_num):
                if i:
                    header = ' '*seq_offset
                else:
                    tag = self[j].tag[:max_tag_len]
                    header = tag+ ' '*(seq_offset-len(tag))
                seq = ' '.join(
                    self.chunks(self[j].seq[s:s+line_len], block_len)
                )
                oh.write(header+seq+'\n')
            oh.write('\n')

    def write_nexus(
        self, oh,
        max_tag_len=10, line_len=60
    ):

        seq_len = len(self[0].seq)
        seq_num = len(self)
        seq_offset = max_tag_len+3

        datatype=self.seq_type().lower()
        if not datatype:
            sys.stderr.write('unknown datatype\n')
            sys.exit(0)
        oh.write('#NEXUS\n\nBegin data;\n')
        oh.write('Dimensions ntax=%s nchar=%s;\n' % (seq_num, seq_len))
        oh.write('Format datatype='+datatype+' interleave=yes gap=-;\n')
        oh.write('Matrix\n')

        for i, s in enumerate(range(0, seq_len, line_len)):
            for j in range(seq_num):
                tag = self[j].tag[:max_tag_len]
                header = tag+' '*(seq_offset-len(tag))
                oh.write(header+self[j].seq[s:s+line_len]+"\n")
            oh.write('\n')
        oh.write('\t;\nend;\n')

    def write_fasta(self, oh, line_len=60):
        for ob in self:
            oh.write(
                ">%s\n%s\n" % (
                    ob.tag,
                    "\n".join(self.chunks(ob.seq, line_len))
                )
            )

    def write(
        self, oh, outfmt='fasta',
        max_tag_len=10, line_len=60, block_len=10, quiet=True
    ):
        outfmt = outfmt.lower()
        if outfmt in ('phylip', 'phy'):
            self.write_phylip(
                oh,
                max_tag_len=max_tag_len,
                line_len=line_len,
                block_len=block_len
            )

            
        elif outfmt == 'tsv':
            for ob in self:
                oh.write(ob.tag+'\t'+ob.seq+'\n')

        elif outfmt == 'csv':
            for ob in self:
                oh.write(ob.tag+','+ob.seq+'\n')

        elif outfmt in ('fasta', 'fas', 'mfa', 'fna', 'fsa', 'fa'):
            self.write_fasta(
                oh,
                line_len=line_len,
            )

        elif outfmt in ('nex', 'nxs', 'nexus'):
            self.write_nexus(
                oh,
                line_len=line_len,
            )

        elif outfmt == 'msf':
            seqs = []
            tags = []
            seq_len = len(self[0].seq)
            seq_offset = max_tag_len + 3
            for o in self:
                tags.append(o.tag[:max_tag_len])
                count = 0
                seq = []
                while 1:
                    frag = o.seq[count:count+block_len]
                    if not frag:
                        break
                    count += block_len
                    seq.append(frag.replace('-', '.'))
                seqs.append(seq)
            block_num = int(line_len/block_len)
            datatype=self.seq_type()
            if datatype == 'Protein':
                datatype = 'P'
            else:
                datatype = 'N'

            oh.write('PileUp\n\n')
            oh.write(' MSF: %s TYPE: %s Check: 0 ..\n\n' % (seq_len, datatype))
            for j in range(len(self)):
                oh.write(
                    ' Name: %s oo Len: %s  Check: 0 Weight: 10.00 ..\n' % (tags[j], len(seqs[j]))
                )
            oh.write('\n//\n\n')
            for i, n in enumerate(range(0, seq_len, line_len)):
                for j in range(len(self)):
                    oh.write(tags[j]+' '*(seq_offset-len(tags[j]))+' '.join(seqs[j][i*block_num:(i+1)*block_num])+'\n')
                oh.write('\n\n')
        
        if not quiet:
            sys.stdout.write('saved at {}\n'.format(oh.name))

    @classmethod
    def parse_string(cls, string, infmt='fasta'):
        o = cls()
        infmt = infmt.lower()
        if infmt in ('fasta', 'fas', 'mfa', 'fna', 'fsa', 'fa'):
            o.parse_fasta(string)
        elif infmt == 'msf':
            o.parse_msf(string)
        elif infmt in ('phylip', 'phy'):
            o.parse_phylip(string)
        return o

    @classmethod
    def parse(cls, fh, infmt='fasta'):
        return cls.parse_string(fh.read(), infmt=infmt)

    @classmethod
    def iterparse(cls, fh):
        o = cls()
        return o.iterparse_fasta(fh)

    @staticmethod
    def chunks(S, n):
        return [S[c:c+n] for c in range(0, len(S), n)]
