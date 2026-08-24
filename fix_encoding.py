from pathlib import Path

p = Path(r"C:\Users\USER\immuno-target-chatbot\app.py")
raw = p.read_bytes()
text = raw.decode("latin-1")
for bad, good in {
    "\x9d": "/",
    "\x97": "-",
    "\x96": "-",
    "\x91": "'",
    "\x92": "'",
    "\x93": '"',
    "\x94": '"',
    "\x85": "...",
    "\x95": "*",
}.items():
    text = text.replace(bad, good)
text = text.replace("Open Targets / Human Protein Atlas / GTEx", "Open Targets / Human Protein Atlas / GTEx")
p.write_text(text, encoding="utf-8")
compile(text, str(p), "exec")
print("app.py encoding and syntax ok")
