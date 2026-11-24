"""`mobilink db`: fetch the three reference databases of the first use case and build the DIAMOND indexes.

CARD (protein homolog model)   the current release from card.mcmaster.ca (its version is written to the manifest)
BacMet2 EXP                    the experimentally confirmed set (753 entries; the packaged class table is keyed to it)
mobileOG-db beatrix-1.6        the archived release on Zenodo (the packaged mobility core is keyed to it)

Each source can also be given as a local file (offline use). The directory gets `mobilink_db.json`, which
`mobilink run` reads: file names, checksums, versions, sources and licences.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import shutil
import subprocess
import tarfile
import urllib.request
from datetime import datetime, timezone
from importlib import resources

from . import __version__

MANIFEST = "mobilink_db.json"
CARD_MEMBER = "protein_fasta_protein_homolog_model.fasta"
SOURCES = {
    "card": dict(
        url="https://card.mcmaster.ca/latest/data",
        licence="CARD licence: free for academic research; commercial use needs a licence from McMaster University "
                "(https://card.mcmaster.ca/about)",
        cite="Alcock et al. 2023, Nucleic Acids Res 51:D690-D699",
    ),
    "bacmet": dict(
        url="http://bacmet.biomedicine.gu.se/download/BacMet2_EXP_database.fasta",
        md5="1043a479c9d94d2a0cfd267147f7691f",
        licence="BacMet 2.0: freely available (http://bacmet.biomedicine.gu.se)",
        cite="Pal et al. 2014, Nucleic Acids Res 42:D737-D743",
    ),
    "mobileog": dict(
        url="https://zenodo.org/api/records/13241605/files/mobileOG-db_beatrix-1.6.All.faa/content",
        md5="5eb498a616e2686b19373f01b0e78eac",
        licence="CC BY 4.0 (Zenodo record 13241605)",
        cite="Brown et al. 2022, Appl Environ Microbiol 88:e00991-22",
    ),
}
UNIVEC_URL = "https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec"


class DbError(RuntimeError):
    """A database could not be fetched, checked or built."""


def data_file(name: str) -> str:
    """Path of a table shipped inside the package (bacmet2_exp_classes.tsv, mobility_core.tsv.gz)."""
    return str(resources.files("mobilink").joinpath("data", name))


def md5sum(path) -> str:
    h = hashlib.md5()  # noqa: S324 (a checksum of a public file, not a security function)
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(url: str, dest: str, log=print) -> None:
    """Download to dest (through dest.part, so an interrupted download is never taken for a file)."""
    req = urllib.request.Request(url, headers={"User-Agent": f"mobilink/{__version__}"})
    log(f"  fetching {url}")
    try:
        with urllib.request.urlopen(req, timeout=120) as r, open(dest + ".part", "wb") as fo:  # noqa: S310
            shutil.copyfileobj(r, fo, 1 << 20)
    except OSError as e:
        raise DbError(f"download failed: {url} ({e}). Download the file by hand and give it to `mobilink db` "
                      "with the matching --card, --bacmet or --mobileog option.") from e
    os.replace(dest + ".part", dest)


def card_from_archive(archive: str, out_dir: str) -> tuple[str, str | None, str]:
    """Protein homolog FASTA, aro_index.tsv (if present) and the CARD version from the CARD data archive."""
    version, index = "unknown", None
    with tarfile.open(archive) as tf:
        names = {os.path.basename(m.name): m for m in tf.getmembers() if m.isfile()}
        if CARD_MEMBER not in names:
            raise DbError(f"{archive}: no {CARD_MEMBER} in the archive")
        for member, target in ((CARD_MEMBER, CARD_MEMBER), ("aro_index.tsv", "aro_index.tsv")):
            if member in names:
                with tf.extractfile(names[member]) as src, open(os.path.join(out_dir, target), "wb") as fo:
                    shutil.copyfileobj(src, fo)
        if "aro_index.tsv" in names:
            index = os.path.join(out_dir, "aro_index.tsv")
        if "card.json" in names:
            with tf.extractfile(names["card.json"]) as src:
                version = str(json.load(src).get("_version", "unknown"))
    return os.path.join(out_dir, CARD_MEMBER), index, version


def mobileog_lookup(fasta: str, out: str) -> int:
    """Lookup table for the mobileOG labels. DIAMOND cuts subject ids at the first space of the FASTA header, and
    many mobileOG-db headers hold a space, so gene, category and evidence are restored from the full header."""
    seen = set()
    op = gzip.open if fasta.endswith(".gz") else open
    with op(fasta, "rt") as fh, gzip.open(out + ".part", "wt") as fo:
        fo.write("sseqid\tgene\tcategory\tevidence\n")
        for ln in fh:
            if ln[0] != ">":
                continue
            h = ln[1:].rstrip("\n")
            key, f = h.split()[0], h.split("|")
            if len(f) < 7:
                raise DbError(f"{fasta}: header without the mobileOG-db fields: {h[:80]}")
            if key in seen:
                raise DbError(f"{fasta}: duplicate DIAMOND key {key}")
            seen.add(key)
            fo.write(f"{key}\t{f[1]}\t{f[3]}\t{f[6]}\n")
    os.replace(out + ".part", out)
    return len(seen)


def count_fasta(path) -> int:
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rb") as fh:
        return sum(1 for ln in fh if ln[:1] == b">")


def tool_version(tool: str, args=("version",)) -> str:
    exe = shutil.which(tool)
    if exe is None:
        raise DbError(f"{tool} is not on the PATH (install it, e.g. `conda install -c bioconda {tool}`)")
    r = subprocess.run([exe, *args], capture_output=True, text=True, check=False)  # noqa: S603
    return (r.stdout or r.stderr).strip().splitlines()[0] if (r.stdout or r.stderr).strip() else "unknown"


def makedb(fasta: str, out_prefix: str, threads: int, diamond: str = "diamond") -> None:
    r = subprocess.run([shutil.which(diamond) or diamond, "makedb", "--in", fasta, "-d", out_prefix, "-p",  # noqa: S603
                        str(threads)], capture_output=True, text=True, check=False)
    if r.returncode != 0 or not os.path.exists(out_prefix + ".dmnd"):
        raise DbError(f"diamond makedb failed for {fasta}:\n{r.stderr[-2000:]}")


def build(out_dir: str, card: str | None = None, bacmet: str | None = None, mobileog: str | None = None,
          threads: int = 4, keep_fasta: bool = False, diamond: str = "diamond", log=print,
          univec: str | bool = True) -> dict:
    """Build the database directory. card: a CARD data archive (.tar.bz2) or the protein homolog FASTA; bacmet and
    mobileog: FASTA files. A source that is None is downloaded. univec: True downloads NCBI UniVec (vector
    sequences, for the check of fixed-length clusters; a failed download is a warning), a path copies that file,
    False leaves it out."""
    os.makedirs(out_dir, exist_ok=True)
    dv = tool_version(diamond)
    info = dict(mobilink=__version__, created=datetime.now(timezone.utc).isoformat(timespec="seconds"), diamond=dv,
                databases={}, tables={})
    work = os.path.join(out_dir, "source")
    os.makedirs(work, exist_ok=True)

    # CARD
    log("CARD (protein homolog model)")
    src, index, version = card, None, "unknown"
    origin = card or SOURCES["card"]["url"]
    if card is None:
        src = os.path.join(work, "card-data.tar.bz2")
        fetch(SOURCES["card"]["url"], src, log)
    if tarfile.is_tarfile(src):
        src, index, version = card_from_archive(src, work)
    makedb(src, os.path.join(out_dir, "card_homolog"), threads, diamond)
    info["databases"]["card"] = dict(dmnd="card_homolog.dmnd", source=origin, version=version,
                                     sequences=count_fasta(src), fasta_md5=md5sum(src),
                                     licence=SOURCES["card"]["licence"], cite=SOURCES["card"]["cite"])
    if index:
        shutil.copyfile(index, os.path.join(out_dir, "aro_index.tsv"))
        info["tables"]["aro_index"] = "aro_index.tsv"

    # BacMet2 EXP
    log("BacMet2 (experimentally confirmed)")
    src, origin = bacmet, bacmet or SOURCES["bacmet"]["url"]
    if bacmet is None:
        src = os.path.join(work, "BacMet2_EXP_database.fasta")
        fetch(SOURCES["bacmet"]["url"], src, log)
    m = md5sum(src)
    note = None
    if m != SOURCES["bacmet"]["md5"]:
        note = "checksum differs from the file the packaged class table was built from"
        log(f"  warning: {note}")
    makedb(src, os.path.join(out_dir, "bacmet2_exp"), threads, diamond)
    shutil.copyfile(data_file("bacmet2_exp_classes.tsv"), os.path.join(out_dir, "bacmet2_exp_classes.tsv"))
    info["databases"]["bacmet"] = dict(dmnd="bacmet2_exp.dmnd", source=origin, version="2.0 EXP",
                                       sequences=count_fasta(src), fasta_md5=m, note=note,
                                       licence=SOURCES["bacmet"]["licence"], cite=SOURCES["bacmet"]["cite"])
    info["tables"]["bacmet_classes"] = "bacmet2_exp_classes.tsv"

    # mobileOG-db
    log("mobileOG-db beatrix-1.6 (355 MB)")
    src, origin = mobileog, mobileog or SOURCES["mobileog"]["url"]
    if mobileog is None:
        src = os.path.join(work, "mobileOG-db_beatrix-1.6.All.faa")
        if not (os.path.exists(src) and md5sum(src) == SOURCES["mobileog"]["md5"]):
            fetch(SOURCES["mobileog"]["url"], src, log)
    m = md5sum(src)
    note = None
    if m != SOURCES["mobileog"]["md5"]:
        if mobileog is None:
            raise DbError("mobileOG-db: the checksum of the download differs from the archived release")
        note = "not the archived beatrix-1.6 file: the packaged mobility core may not match its entries"
        log(f"  warning: {note}")
    n = mobileog_lookup(src, os.path.join(out_dir, "mobileog_all.lookup.tsv.gz"))
    makedb(src, os.path.join(out_dir, "mobileog_all"), threads, diamond)
    shutil.copyfile(data_file("mobility_core.tsv.gz"), os.path.join(out_dir, "mobility_core.tsv.gz"))
    info["databases"]["mobileog"] = dict(dmnd="mobileog_all.dmnd", source=origin, version="beatrix-1.6",
                                         sequences=n, fasta_md5=m, note=note,
                                         licence=SOURCES["mobileog"]["licence"], cite=SOURCES["mobileog"]["cite"])
    info["tables"].update(mobileog_lookup="mobileog_all.lookup.tsv.gz", mobility_core="mobility_core.tsv.gz")

    if univec:
        log("NCBI UniVec (vector sequences)")
        dest = os.path.join(out_dir, "UniVec.fa")
        try:
            if univec is True:
                fetch(UNIVEC_URL, dest, log)
            else:
                shutil.copyfile(univec, dest)
            info["tables"]["univec"] = "UniVec.fa"
            info["univec"] = dict(source=UNIVEC_URL if univec is True else univec, sequences=count_fasta(dest),
                                  md5=md5sum(dest), licence="public domain (NCBI)")
        except DbError as e:
            log(f"  warning: UniVec not installed ({e}); fixed-length clusters will be flagged but not checked")

    if not keep_fasta:
        shutil.rmtree(work, ignore_errors=True)
    with open(os.path.join(out_dir, MANIFEST), "w") as fh:
        json.dump(info, fh, indent=2)
    return info


def load(db_dir: str) -> dict:
    """The manifest of a database directory, with absolute paths under `paths`; every file is checked to exist."""
    path = os.path.join(db_dir, MANIFEST)
    if not os.path.exists(path):
        raise DbError(f"{db_dir} is not a mobilink database directory (no {MANIFEST}); build one with `mobilink db`")
    with open(path) as fh:
        info = json.load(fh)
    paths = {k: os.path.join(db_dir, v["dmnd"]) for k, v in info["databases"].items()}
    paths.update({k: os.path.join(db_dir, v) for k, v in info["tables"].items()})
    missing = [p for p in paths.values() if not os.path.exists(p)]
    if missing:
        raise DbError(f"{db_dir}: files named in {MANIFEST} are missing: {', '.join(missing)}")
    info["paths"] = paths
    return info
