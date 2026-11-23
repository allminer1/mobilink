"""P4 step 4.4 (kit rule, record part): what the archive records say about the library kit and chemistry of each
project of the sample. One selected run per project and class is looked up at NCBI (efetch, db=sra): the design
description, the library construction protocol, the instrument model and any run attribute that names a kit, flow
cell or basecaller. Metadata only; no reads. Output: public/screen/p4_record_kits.tsv (for reading; the kit table
public/screen/p4_project_kits.tsv is entered by hand from this, the read-end probe and the publications).
usage: python -I P4_02b_record_kits.py
"""
import re
import time
import urllib.request

import pandas as pd

L = "/mnt/rmrfc/workflow/iskander/linkstat"
KEY = re.compile(r"(SQK-[A-Z0-9.\-]+|FLO-[A-Z0-9]+|R ?9\.?\d(?:\.\d)?|R ?10\.?\d(?:\.\d)?|LSK ?\d{3}|RBK ?\d{3}|NBD ?\d{3}|RAD ?\d{3}|RPB ?\d{3}|PCR|"
                 r"ligation|rapid|native barcod\w*|guppy[ v\d.]*|dorado[ v\d.]*|albacore[ v\d.]*|bonito[ v\d.]*|sup\b|hac\b|kit ?\d{1,3}\b|MDA|amplif\w+)", re.I)


def get(url, tries=3):
    for k in range(tries):
        try:
            return urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "metadata-screen"}), timeout=60).read().decode(errors="replace")
        except Exception as e:  # noqa: BLE001
            err = str(e)
            time.sleep(4 * (k + 1))
    return f"ERROR {err}"


def tag(x, name):
    m = re.search(rf"<{name}[^>]*>(.*?)</{name}>", x, re.S)
    return re.sub(r"\s+", " ", m.group(1)).strip() if m else ""


def main():
    s = pd.read_csv(f"{L}/public/screen/p4_selected.tsv", sep="\t", dtype=str)
    rows = []
    for (p, c), g in s.groupby(["study_accession", "class"]):
        run = g.run_accession.iloc[0]
        x = get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id={run}&rettype=xml")
        time.sleep(0.5)
        design, proto = tag(x, "DESIGN_DESCRIPTION"), tag(x, "LIBRARY_CONSTRUCTION_PROTOCOL")
        attrs = "; ".join(f"{a}={b}" for a, b in re.findall(r"<TAG>(.*?)</TAG>\s*<VALUE>(.*?)</VALUE>", x) if KEY.search(a + " " + b))[:400]
        text = " | ".join(t for t in (design, proto, attrs, tag(x, "TITLE")) if t)
        rows.append(dict(project=p, env=c, run=run, runs_selected=len(g), instrument=tag(x, "INSTRUMENT_MODEL"), first_public=g.first_public.min(),
                         hits="; ".join(sorted({m.group(0).strip() for m in KEY.finditer(text)}, key=str.lower)),
                         design=design[:500], protocol=proto[:500], attributes=attrs, error=x[:80] if x.startswith("ERROR") else ""))
    t = pd.DataFrame(rows)
    t.to_csv(f"{L}/public/screen/p4_record_kits.tsv", sep="\t", index=False)
    print(len(t), "project-class rows;", int((t.hits != "").sum()), "with a kit, chemistry or basecaller word;", int((t.error != "").sum()), "errors")


if __name__ == "__main__":
    main()
