import os, random, subprocess, sys, tempfile
AO = os.path.abspath(sys.argv[1])
os.chdir(tempfile.mkdtemp())
random.seed(1)
rows = []
for i in range(100000):   # abundant dispersed family: 100k copies over 2 Gb
    c = "ctg%02d" % random.randrange(20); s = random.randrange(1, 100_000_000)
    rows.append(["RHIN", str(random.randint(100, 250)), c, str(s), str(s + 180), "+"])
for i in range(20):       # one tandem array, near-identical units, top scores
    s = 5_000_000 + i * 3000
    rows.append(["RHIN", str(400 - i), "ctg03", str(s), str(s + 180), "+"])
rows.sort(key=lambda r: -int(r[1]))
open("t_sorted.tsv", "w").write("".join("\t".join(r) + "\n" for r in rows))
random.shuffle(rows)
open("t_shuf.tsv", "w").write("".join("\t".join(r) + "\n" for r in rows))
top = [l.split("\t") for l in subprocess.run([sys.executable, AO, "t_sorted.tsv"], capture_output=True, text=True).stdout.splitlines()][:100]
rnd = [l.split("\t") for l in subprocess.run([sys.executable, AO, "t_shuf.tsv", "--mark-only", "--limit", "100"], capture_output=True, text=True).stdout.splitlines()][:100]
arr = lambda r: r[7] == "array"
unit = lambda r: r[2] == "ctg03" and 5_000_000 <= int(r[3]) < 5_060_000
print("top100: rows %d, marked array %d, array units on plate %d (want 1), first row is an array unit: %s"
      % (len(top), sum(map(arr, top)), sum(map(unit, top)), unit(top[0])))
print("top100: dispersed rows marked array %d (want 0)" % sum(1 for r in top if arr(r) and not unit(r)))
print("rand100: marked array %d (want ~0; array units drawn %d)" % (sum(map(arr, rnd)), sum(map(unit, rnd))))
