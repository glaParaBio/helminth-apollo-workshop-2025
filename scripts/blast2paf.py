#!/usr/bin/env python3

import sys

if len(sys.argv) > 1 and ('-h' in sys.argv or '--help' in sys.argv):
    print("""\
Convert tabular blast output to paf format. See the code itself for expected
columns from blast
USAGE
cat blast.out | blast2paf.py > blast.paf""")
    sys.exit(0)

for line in sys.stdin:
    if line.startswith('#'):
        continue
    # These are the columns we expect in blast input
    [qaccver, qlen, qstart, qend, sstrand, saccver, slen, sstart, send, nident] = line.strip().split('\t')
    qlen = int(qlen)
    qstart = int(qstart)
    qend = int(qend)
    slen = int(slen)
    sstart = int(sstart)
    send = int(send)
    nident = int(nident)
    
    if ((qstart < qend) and (sstart < send)) or ((qstart > qend) and (sstart > send)):
        sstrand = '+'
    else:
        sstrand = '-'

    if qstart > qend:
        _qstart = qend
        _qend = qstart
        qstart = _qstart
        qend = _qend
    if sstart > send:
        _sstart = send
        _send = sstart
        sstart = _sstart
        send = _send
    qstart -= 1
    sstart -= 1
    aln_length = qend - qstart
    mapq = 255
    out = [qaccver, qlen, qstart, qend, sstrand, saccver, slen, sstart, send, nident, aln_length, mapq]
    out = '\t'.join([str(x) for x in out])
    print(out)
