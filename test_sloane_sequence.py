import re
import gzip
from time import time

import sloane_sequence

IMPLEMENTED_SEQUENCES = [i for i in dir(
    sloane_sequence) if re.match("A[0-9]{6}", i)]


OEIS_DB_PATH = "./oeis_data/stripped.gz" # downloaded from https://oeis.org/stripped.gz
OEIS_NAMES_PATH = "./oeis_data/names.gz" # downloaded from https://oeis.org/names.gz

def test_seq_terms():
    for line in gzip.open(OEIS_DB_PATH):
        if line[0] == b"#"[0]:
            continue
        else:
            seq_name = line.split()[0].decode()
            if seq_name in IMPLEMENTED_SEQUENCES:
                actual_seq_terms = [int(i) for i in line.split()[
                    1].split(b",")[1:-1]]
                start = time()
                seq_obj = getattr(sloane_sequence, seq_name)()
                implemented_seq_terms = seq_obj._eval_up_to_n(
                    len(actual_seq_terms))
                if actual_seq_terms != implemented_seq_terms:
                    print(seq_name, implemented_seq_terms, actual_seq_terms, len(
                        implemented_seq_terms), len(actual_seq_terms))
                    for i, (implemented, actual) in enumerate(zip(implemented_seq_terms, actual_seq_terms)):
                        if implemented != actual:
                            raise ValueError(
                                f"{seq_name}({seq_obj.offset+i}) != {implemented}, should be {actual}")
                else:
                    end = time()
                    print(f"[✓] {seq_name}, time={end-start}s")
    print("[+] All generated terms correct!")


def test_fields():
    for line in gzip.open(OEIS_NAMES_PATH):
        if line[0] == b"#"[0]:
            continue
        else:
            seq_name = line[:7].decode()
            seq_description = line[8:].strip().decode()
            if seq_name in IMPLEMENTED_SEQUENCES:
                seq_obj = getattr(sloane_sequence, seq_name)()
                if seq_obj.seq_number != int(seq_name[1:]):
                    raise ValueError(f"wrong sequence number: {seq_name}.seq_number is {seq_obj.seq_number}, should be {int(seq_name[1:])}")
                elif seq_description != seq_obj.description:
                    raise ValueError(f"wrong description: {seq_name}.description is {repr(seq_obj.description)}, should be {repr(seq_description)}")
                else:
                    print(f"[✓] {seq_name}")
    print("[+] All sequences have correct information in fields")

test_seq_terms()
test_fields()