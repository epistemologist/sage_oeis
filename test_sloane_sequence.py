import re
import bz2
from tqdm import tqdm
from time import time

import sloane_sequence

OEIS_DB_PATH = "/usr/share/sloane/sloane-oeis.bz2"
OEIS_DB = [line for line in bz2.open(OEIS_DB_PATH).read().splitlines() if line[0] != ord("#")]

print("[+] Reading sequence data...")
SEQUENCE_DATA = {
    line.split()[0].decode() : [int(i) for i in line.split()[1].split(b",")[1:-1]]
    for line in tqdm(OEIS_DB)
}

IMPLEMENTED_SEQUENCES = [ i for i in dir(sloane_sequence) if re.match("A[0-9]{6}", i)]

for seq in sorted(IMPLEMENTED_SEQUENCES):
    start = time()
    seq_obj = getattr(sloane_sequence, seq)()
    actual_sequence = SEQUENCE_DATA[seq]
    implemented_sequence = seq_obj._eval_up_to_n(len(actual_sequence))
    if actual_sequence != implemented_sequence:
        print(seq, implemented_sequence, actual_sequence, len(implemented_sequence), len(actual_sequence))
        for i, (implemented, actual) in enumerate(zip(implemented_sequence, actual_sequence)):
            if implemented != actual:
                raise ValueError(f"{seq}({seq_obj.offset + i}) != {implemented}, should be {actual}")
    else:
        end = time()
        print(f"[✓] {seq}, time={end-start}s")
    # print(seq, actual_sequence, implemented_sequence)