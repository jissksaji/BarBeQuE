import hyperscan

sequence = "TTTTATCGGAAAAAACCCCGGTACCTTTT"

forward_primer = "ATCGGA"
reverse_primer_rc = "GGTACC"

patterns = [
    forward_primer.encode(),
    reverse_primer_rc.encode()
]

db = hyperscan.Database()
db.compile(
    expressions=patterns,
    ids=[0, 1],
    elements=2
)

matches = []

def on_match(id, from_, to, flags, context):
    name = "F" if id == 0 else "R_rc"
    matches.append((name, from_, to))
    return 0

db.scan(sequence.encode(), match_event_handler=on_match)

print(matches)