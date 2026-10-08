from pathlib import Path
import hashlib

base = Path(__file__).resolve().parent
count = 0
for line in (base / 'SHA256SUMS.txt').read_text(encoding='utf-8').splitlines():
    expected, name = line.split('  ', 1)
    actual = hashlib.sha256((base / name).read_bytes()).hexdigest()
    if actual != expected:
        raise SystemExit('HASH MISMATCH: ' + name)
    count += 1
print('All ' + str(count) + ' archived file hashes match.')
