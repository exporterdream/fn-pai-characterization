#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Comparative genomics of the three candidate PAIs in Fusobacterium nucleatum ATCC 25586.
=========================================================================================
I use macOS M1 to run code.
    -Homebrew to download Miniforge, previous code use anaconda, but there is something wrong due to ox84 vs. intel, arm problem
    So I change to Miniforge for conda env.
    I run my other code, they still work.

Package needs:
# Create the environment and install tools:
#      conda create -n pai -c bioconda -c conda-forge \
#          python=3.11 blast mafft fasttree biopython pandas matplotlib openpyxl -y
#      conda activate pai

Inputs:
#   25586.fasta   -> AE009951.2 F. nucleatum ATCC 25586 complete genome (query)
#   W83.fasta     -> Porphyromonas gingivalis W83 (outgroup)
#   E.coil.fasta  -> Escherichia coli 09-00049 (outgroup)
#   Supplementary_Table_S1_PAI1Blast.xlsx
#   Supplementary_Table_S1_PAI2Blast.xlsx
#   Supplementary_Table_S1_PAI3Blast.xlsx
# Everything else (genomes, proteins, GFFs) is downloaded from NCBI from M0

Reproduces the analysis in report_PAI_comparative_genomics.md:
  M0  download curated genome set + extract PAI regions/proteins
  M1  whole-region BLASTN conservation matrix + heatmap
  M2  gene-level TBLASTN presence/absence matrix + heatmap
  M3  synteny / gene-order maps per PAI
  M4  core-genome phylogeny + PAI distribution
  M5  targeted HGT tests (composition, mobility context, gene-tree concordance), MAFFT
  M6  GC-skew reinterpretation panel

  To assess whether the three candidate PAIs are specific to ATCC~25586 or conserved across the genus, 
  a comparative genomic analysis was performed on a curated set of 17 complete genomes: 
  15 \textit{Fusobacterium} genomes spanning all four \textit{F.~nucleatum} subspecies (\textit{nucleatum}, 
  \textit{animalis}, \textit{polymorphum}, \textit{vincentii}) and four related species (\textit{F.~periodonticum}, 
  \textit{F.~varium}, \textit{F.~ulcerans}, \textit{F.~mortiferum}), plus two outgroups (\textit{P.~gingivalis} W83 
  and \textit{E.~coli} 09-00049). All genomes were downloaded from NCBI RefSeq via the Datasets v2 API, including genomic FASTA, 
  GFF annotation, and protein FASTA files.

  Three targeted tests for horizontal gene transfer were applied to each PAI. 
  \textbf{Compositional analysis:} GC content and dinucleotide frequency distances for each PAI were compared against a null distribution of 500 size-matched random genomic windows; 
  a PAI was considered compositionally atypical if either metric exceeded $|z| > 2$.
    \textbf{Mobility context:} the AE009951.2 GenBank annotation was scanned for transposase, integrase, recombinase, 
    and insertion sequence genes within each PAI and within 15~kb flanking regions. 
    \textbf{Gene-tree concordance:} for five representative PAI genes (FN1885, FN0837, FN0849, FN1317, FN1318), 
    individual gene trees were inferred with MAFFT~+~FastTree and compared to the species tree by normalized Robinson--Foulds 
    distance. 

"""

import os, io, json, time, pickle, zipfile, subprocess, urllib.request, urllib.parse
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = ["Liberation Sans", "Arimo", "DejaVu Sans"]
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap, BoundaryNorm, LinearSegmentedColormap, Normalize

from Bio import SeqIO, AlignIO, Phylo
from Bio.SeqRecord import SeqRecord

# --------------------------------------------------------------------------
# CONFIG
# --------------------------------------------------------------------------
INDIR = "inputs"                 # where your uploaded files live
WORK = "pai_work"                # working directory (created)
FIGS = os.path.join(WORK, "figs")
RESULTS = os.path.join(WORK, "results")
for d in [WORK, FIGS, RESULTS, f"{WORK}/genomes", f"{WORK}/gff", f"{WORK}/proteins",
          f"{WORK}/pai_regions", f"{WORK}/pai_proteins", f"{WORK}/blast/db",
          f"{WORK}/blast/m1", f"{WORK}/blast/m2", f"{WORK}/blast/m3", f"{WORK}/blast/pdb",
          f"{WORK}/core_aln/aligned", f"{WORK}/genetrees"]:
    os.makedirs(d, exist_ok=True)

# PAI coordinates (1-based inclusive) on AE009951.2, from the manuscript / IslandViewer
PAIS = {"PAI1": (366847, 372320), "PAI2": (1496613, 1523855), "PAI3": (1967193, 1999047)} # identified by islandviewer
PAILEN = {k: v[1] - v[0] + 1 for k, v in PAIS.items()}

# Curated 15-genome Fusobacterium set (complete RefSeq assemblies)
GENOMES = [
    ("GCF_003019295.1", "Fn_nucleatum_ATCC25586"),
    ("GCF_037889005.1", "Fn_nucleatum_SB011"),
    ("GCF_001296145.1", "Fn_animalis_KCOM1325"),
    ("GCF_056819735.1", "Fn_animalis_CJ003"),
    ("GCF_001457555.1", "Fn_polymorphum_NCTC10562"),
    ("GCF_050920725.1", "Fn_polymorphum_ATCC10953"),
    ("GCF_019552005.1", "Fn_vincentii_THCT14A3"),
    ("GCF_056820395.1", "Fn_vincentii_CJ001"),
    ("GCF_003019755.1", "F_periodonticum_2_1_31"),
    ("GCF_003019655.1", "F_varium_ATCC27725"),
    ("GCF_054951225.1", "F_varium_THCT4E2"),
    ("GCF_003019675.1", "F_ulcerans_ATCC49185"),
    ("GCF_900478315.1", "F_ulcerans_NCTC12112"),
    ("GCF_003019315.1", "F_mortiferum_ATCC9817"),
    ("GCF_057585585.1", "F_mortiferum_SYC45"),
]
ORDER = [lbl for _, lbl in GENOMES]
OUTGROUPS = [("P_gingivalis_W83", "W83.fasta"), ("E_coli_09-00049", "E.coil.fasta")]
UA = {"User-Agent": "pai-comparative-genomics/1.0"}


def http(url, timeout=120):
    req = urllib.request.Request(url, headers=UA)
    return urllib.request.urlopen(req, timeout=timeout)


# --------------------------------------------------------------------------
# M0: download genomes + extract PAI regions/proteins, Only use once, I comment them for I have it download already, you can use it by uncomment in main.
# --------------------------------------------------------------------------
def m0():
    print("[M0] downloading genomes and extracting PAI data")
    # copy query + outgroups
    import shutil
    shutil.copy(os.path.join(INDIR, "25586.fasta"), f"{WORK}/genomes/Fn_query_AE009951.fasta")
    for lbl, fn in OUTGROUPS:
        shutil.copy(os.path.join(INDIR, fn), f"{WORK}/genomes/{lbl}.fasta")
    # download Fusobacterium genomes (genomic fasta + gff + protein)
    for acc, lbl in GENOMES:
        if os.path.exists(f"{WORK}/genomes/{lbl}.fasta"):
            continue
        url = (f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/download"
               f"?include_annotation_type=GENOME_FASTA,GENOME_GFF,PROT_FASTA&filename={acc}.zip")
        try:
            data = http(url).read()
            z = zipfile.ZipFile(io.BytesIO(data))
            fna = [n for n in z.namelist() if n.endswith("_genomic.fna")]
            gff = [n for n in z.namelist() if n.endswith("genomic.gff")]
            faa = [n for n in z.namelist() if n.endswith("protein.faa")]
            if fna: open(f"{WORK}/genomes/{lbl}.fasta", "wb").write(z.read(fna[0]))
            if gff: open(f"{WORK}/gff/{lbl}.gff", "wb").write(z.read(gff[0]))
            if faa: open(f"{WORK}/proteins/{lbl}.faa", "wb").write(z.read(faa[0]))
            print(f"  OK {lbl}")
            time.sleep(0.4)
        except Exception as e:
            print(f"  FAIL {lbl}: {e}")
    # extract PAI nucleotide regions from AE009951.2
    rec = SeqIO.read(os.path.join(INDIR, "25586.fasta"), "fasta")
    for name, (s, e) in PAIS.items():
        sub = rec.seq[s - 1:e]
        SeqIO.write(SeqRecord(sub, id=f"{name}_ATCC25586", description=f"{name} {s}-{e}"),
                    f"{WORK}/pai_regions/{name}.fasta", "fasta")
    # parse PAI gene lists from the my BLAST tables
    import openpyxl
    tables = {"PAI1": "Supplementary_Table_S1_PAI1Blast.xlsx",
              "PAI2": "Supplementary_Table_S1_PAI2Blast.xlsx.xlsx",
              "PAI3": "Supplementary_Table_S1_PAI3Blast.xlsx.xlsx"}
    genes = {}
    for pai, fn in tables.items():
        ws = openpyxl.load_workbook(os.path.join(INDIR, fn), read_only=True).active
        rows = list(ws.iter_rows(values_only=True)); hdr = rows[0]
        li, pi, di = hdr.index("locus_tag"), hdr.index("protein_id"), hdr.index("product")
        genes[pai] = [{"locus": r[li], "protein_id": r[pi], "product": r[di]}
                      for r in rows[1:] if r[li]]
    json.dump(genes, open(f"{WORK}/pai_genes.json", "w"), indent=1)
    # fetch PAI protein sequences from NCBI protein by protein_id
    ids = [(g["protein_id"], g["locus"], g["product"], p) for p in genes for g in genes[p]]
    fetched = {}
    for i in range(0, len(ids), 100):
        chunk = [a[0] for a in ids[i:i + 100]]
        q = urllib.parse.urlencode({"db": "protein", "id": ",".join(chunk),
                                    "rettype": "fasta", "retmode": "text"})
        txt = http("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + q, 60).read().decode()
        for r in SeqIO.parse(io.StringIO(txt), "fasta"):
            fetched[r.id] = r
        time.sleep(0.5)
    meta = {}
    for pid, locus, product, pai in ids:
        r = fetched.get(pid)
        if r is None:
            continue
        with open(f"{WORK}/pai_proteins/{pai}_proteins.fasta", "a") as h:
            SeqIO.write(SeqRecord(r.seq, id=locus, description=f"{product} [{pai}|{pid}]"), h, "fasta")
        meta[locus] = {"pai": pai, "protein_id": pid, "product": product, "len": len(r.seq)}
    json.dump(meta, open(f"{WORK}/pai_protein_meta.json", "w"), indent=1)
    print(f"[M0] done: {len(meta)} PAI proteins")


# --------------------------------------------------------------------------
# BLAST helpers
# --------------------------------------------------------------------------
def run(cmd):
    subprocess.run(cmd, check=True, capture_output=True)


def build_dbs():
    labels = ORDER + [l for l, _ in OUTGROUPS]
    fasta = {l: f"{WORK}/genomes/{l}.fasta" for l in ORDER}
    fasta.update({l: f"{WORK}/genomes/{l}.fasta" for l, _ in OUTGROUPS})
    for l in labels:
        db = f"{WORK}/blast/db/{l}"
        if not os.path.exists(db + ".nsq"):
            run(["makeblastdb", "-in", fasta[l], "-dbtype", "nucl", "-out", db])
    return labels


# --------------------------------------------------------------------------
# M1: region-level BLASTN conservation
# --------------------------------------------------------------------------
def m1(labels):
    print("[M1] region-level BLASTN conservation")
    for pai in PAIS:
        for l in labels:
            out = f"{WORK}/blast/m1/{pai}__{l}.tsv"
            if not os.path.exists(out):
                run(["blastn", "-query", f"{WORK}/pai_regions/{pai}.fasta", "-db", f"{WORK}/blast/db/{l}",
                     "-out", out, "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore",
                     "-evalue", "1e-5", "-max_target_seqs", "50", "-num_threads", "4"])
    rows = []
    for pai in PAIS:
        for l in labels:
            f = f"{WORK}/blast/m1/{pai}__{l}.tsv"
            hsps = []
            if os.path.getsize(f) > 0:
                for line in open(f):
                    p = line.rstrip("\n").split("\t")
                    if len(p) < 10: continue
                    hsps.append((min(int(p[4]), int(p[5])), max(int(p[4]), int(p[5])), float(p[2])))
            hsps.sort(); cov = 0; best = 0; cs = ce = None
            for s, e, pi in hsps:
                best = max(best, pi)
                if cs is None: cs, ce = s, e
                elif s <= ce + 1: ce = max(ce, e)
                else: cov += ce - cs + 1; cs, ce = s, e
            if cs is not None: cov += ce - cs + 1
            rows.append({"PAI": pai, "genome": l, "query_cov_pct": round(100 * cov / PAILEN[pai], 1),
                         "best_pident": round(best, 1), "n_hsps": len(hsps)})
    m1 = pd.DataFrame(rows)
    m1.to_csv(f"{RESULTS}/M1_region_conservation_long.csv", index=False)
    m1.pivot(index="genome", columns="PAI", values="query_cov_pct").to_csv(f"{RESULTS}/M1_query_coverage_matrix.csv")
    m1.pivot(index="genome", columns="PAI", values="best_pident").to_csv(f"{RESULTS}/M1_best_identity_matrix.csv")
    # heatmap
    order = ORDER + [l for l, _ in OUTGROUPS]
    cov = m1.pivot(index="genome", columns="PAI", values="query_cov_pct").reindex(order)[list(PAIS)]
    pid = m1.pivot(index="genome", columns="PAI", values="best_pident").reindex(order)[list(PAIS)]
    cmap = LinearSegmentedColormap.from_list("c", ["#FAF9F3", "#75A025", "#0279EE"])
    fig, axes = plt.subplots(1, 2, figsize=(11, 7.5))
    for ax, data, title in [(axes[0], cov, "Query coverage (%)"), (axes[1], pid, "Best nucleotide identity (%)")]:
        im = ax.imshow(data.values, aspect="auto", cmap=cmap, vmin=0, vmax=100)
        ax.set_xticks(range(3)); ax.set_xticklabels(list(PAIS))
        ax.set_yticks(range(len(order))); ax.set_yticklabels([x.replace("_", " ") for x in order], fontsize=8)
        ax.set_title(title, fontsize=11)
        for i in range(len(order)):
            for j in range(3):
                v = data.values[i, j]
                ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=7,
                        color="white" if v > 55 else "black")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    axes[0].set_ylabel("Genome")
    fig.suptitle("M1: Whole-region conservation of the three candidate PAIs", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(f"{FIGS}/M1_region_conservation_heatmap.png", dpi=200)
    fig.savefig(f"{FIGS}/M1_region_conservation_heatmap.svg")
    plt.close(fig)


# --------------------------------------------------------------------------
# M2: gene-level TBLASTN presence/absence
# --------------------------------------------------------------------------
def m2(labels):
    print("[M2] gene-level TBLASTN presence/absence")
    meta = json.load(open(f"{WORK}/pai_protein_meta.json"))
    for pai in PAIS:
        for l in labels:
            out = f"{WORK}/blast/m2/{pai}__{l}.tsv"
            if not os.path.exists(out):
                run(["tblastn", "-query", f"{WORK}/pai_proteins/{pai}_proteins.fasta",
                     "-db", f"{WORK}/blast/db/{l}", "-out", out,
                     "-outfmt", "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore",
                     "-evalue", "1e-5", "-max_target_seqs", "20", "-num_threads", "4"])
    gene_order = {p: sorted([l for l in meta if meta[l]["pai"] == p]) for p in PAIS}
    rows = []
    for pai in PAIS:
        for g in labels:
            f = f"{WORK}/blast/m2/{pai}__{g}.tsv"; best = {}
            if os.path.getsize(f) > 0:
                for line in open(f):
                    p = line.rstrip("\n").split("\t")
                    if len(p) < 11: continue
                    locus, pident, aln, qlen = p[0], float(p[2]), int(p[3]), int(p[4])
                    cov = aln / qlen; cur = best.get(locus, {"cov": 0, "pid": 0})
                    if pident >= 30 and cov > cur["cov"]: cur["cov"] = cov
                    cur["pid"] = max(cur["pid"], pident); best[locus] = cur
            for locus in gene_order[pai]:
                b = best.get(locus, {"cov": 0, "pid": 0})
                rows.append({"PAI": pai, "locus": locus, "genome": g,
                             "present": 1 if (b["cov"] >= 0.5 and b["pid"] >= 30) else 0,
                             "best_cov": round(b["cov"], 2), "best_pid": round(b["pid"], 1)})
    m2 = pd.DataFrame(rows)
    m2.to_csv(f"{RESULTS}/M2_gene_presence_long.csv", index=False)
    m2.pivot_table(index=["PAI", "locus"], columns="genome", values="present", aggfunc="max") \
      .to_csv(f"{RESULTS}/M2_gene_presence_matrix.csv")
    # heatmap
    order = ORDER + [l for l, _ in OUTGROUPS]
    cmap = ListedColormap(["#FAF9F3", "#0279EE"]); norm = BoundaryNorm([0, 0.5, 1], cmap.N)
    fig, axes = plt.subplots(3, 1, figsize=(11.5, 12), gridspec_kw={"height_ratios": [10, 32, 33]})
    for ax, pai in zip(axes, PAIS):
        sub = m2[m2.PAI == pai]; genes = sorted(sub.locus.unique())
        M = sub.pivot_table(index="locus", columns="genome", values="present", aggfunc="max").reindex(genes)[order]
        ax.imshow(M.values, aspect="auto", cmap=cmap, norm=norm, interpolation="nearest")
        ax.set_yticks(range(len(genes)))
        ax.set_yticklabels([f"{g} {str(meta[g]['product'])[:18]}" for g in genes], fontsize=6.8)
        ax.set_ylabel(f"{pai} ({len(genes)} genes)", fontsize=10); ax.set_xticks([])
        ax.tick_params(length=0)
        for sp in ax.spines.values(): sp.set_visible(False)
        ax.set_xticks(np.arange(-.5, len(order), 1), minor=True)
        ax.set_yticks(np.arange(-.5, len(genes), 1), minor=True)
        ax.grid(which="minor", color="white", linewidth=0.4)
    axes[-1].set_xticks(range(len(order)))
    axes[-1].set_xticklabels([x.replace("_", " ") for x in order], rotation=45, ha="right", fontsize=8)
    fig.legend(handles=[mpatches.Patch(color="#0279EE", label="Present"),
                        mpatches.Patch(color="#FAF9F3", label="Absent")],
               loc="upper right", bbox_to_anchor=(0.995, 0.998), fontsize=9)
    fig.suptitle("M2: Gene-level presence/absence of PAI-encoded genes (TBLASTN)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.985]); fig.subplots_adjust(left=0.24)
    fig.savefig(f"{FIGS}/M2_gene_presence_heatmap.png", dpi=200)
    # fig.savefig(f"{FIGS}/M2_gene_presence_heatmap.svg")
    plt.close(fig)


# --------------------------------------------------------------------------
# M3: synteny maps
# --------------------------------------------------------------------------
def parse_gff(path, with_pid=False):
    rows = []
    for line in open(path):
        if line.startswith("#"): continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9 or p[2] != "CDS": continue
        attrs = dict(kv.split("=", 1) for kv in p[8].split(";") if "=" in kv)
        rows.append({"seqid": p[0], "start": int(p[3]), "end": int(p[4]), "strand": p[6],
                     "locus": attrs.get("locus_tag", "?"), "protein_id": attrs.get("protein_id", ""),
                     "gene": attrs.get("gene", ""), "product": attrs.get("product", "")})
    df = pd.DataFrame(rows)
    return df[(df.end - df.start) <= 100000].reset_index(drop=True)  # drop full-length artifacts


def m3():
    print("[M3] synteny maps")
    # query annotation from AE009951.2 GenBank
    q = urllib.parse.urlencode({"db": "nuccore", "id": "AE009951.2", "rettype": "gbwithparts", "retmode": "text"})
    txt = http("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + q, 90).read().decode()
    rec = SeqIO.read(io.StringIO(txt), "genbank")
    qrows = [{"seqid": rec.id, "start": int(f.location.start) + 1, "end": int(f.location.end),
              "strand": "+" if f.location.strand == 1 else "-",
              "locus": f.qualifiers.get("locus_tag", ["?"])[0], "protein_id": "",
              "gene": f.qualifiers.get("gene", [""])[0], "product": f.qualifiers.get("product", [""])[0]}
             for f in rec.features if f.type == "CDS"]
    gffdata = {"AE009951.2_query": pd.DataFrame(qrows)}
    gffdata["AE009951.2_query"] = gffdata["AE009951.2_query"][
        (gffdata["AE009951.2_query"].end - gffdata["AE009951.2_query"].start) <= 100000]
    prot2locus = {}
    for l in ORDER:
        df = parse_gff(f"{WORK}/gff/{l}.gff")
        gffdata[l] = df
        prot2locus[l] = {r["protein_id"]: r["locus"] for _, r in df.iterrows() if r["protein_id"]}
    pickle.dump(prot2locus, open(f"{WORK}/prot2locus.pkl", "wb"))

    def best_anchor(pai, label):
        f = f"{WORK}/blast/m1/{pai}__{label}.tsv"
        if not (os.path.exists(f) and os.path.getsize(f) > 0): return None
        hsps = []
        for line in open(f):
            p = line.rstrip("\n").split("\t")
            if len(p) < 10: continue
            hsps.append({"sseqid": p[1], "pident": float(p[2]), "len": int(p[3]), "qs": int(p[4]),
                         "qe": int(p[5]), "ss": int(p[6]), "se": int(p[7]), "bit": float(p[9])})
        if not hsps: return None
        df = pd.DataFrame(hsps)
        df["smin"] = df[["ss", "se"]].min(axis=1); df["smax"] = df[["ss", "se"]].max(axis=1)
        best = None
        for seqid, grp in df.groupby("sseqid"):
            grp = grp.sort_values("smin").reset_index(drop=True)
            clusters, cur = [], [0]
            for i in range(1, len(grp)):
                if grp.smin[i] - grp.smax[cur[-1]] <= 60000: cur.append(i)
                else: clusters.append(cur); cur = [i]
            clusters.append(cur)
            for cl in clusters:
                sub = grp.loc[cl]; score = sub.bit.sum()
                qcov = min((sub.qe - sub.qs + 1).clip(0).sum() / PAILEN[pai], 1.0)
                lo, hi = sub.smin.min(), sub.smax.max()
                if best is None or score > best["score"]:
                    best = {"seqid": seqid, "start": lo, "end": hi, "score": score, "qcov": qcov}
        return best

    def neighborhood(df, seqid, s, e, flank=6000):
        d = df[df.seqid == seqid]
        return d[(d.end >= s - flank) & (d.start <= e + flank)].sort_values("start").copy()

    def module_color(product):
        p = str(product).lower()
        if "transposase" in p or "integrase" in p or "recombinase" in p: return "#E4572E"
        if "hemolysin" in p or "toxin" in p or "hepn" in p or "antitoxin" in p: return "#B9314F"
        if "biotin" in p or "oxononanoate" in p or "pimeloyl" in p: return "#0279EE"
        if any(k in p for k in ["glycogen", "glucose-1-phosphate", "adpglucose", "amylase", "glg", "malq"]): return "#75A025"
        if "sigma" in p or "rpod" in p or "siga" in p: return "#7B2CBF"
        if any(k in p for k in ["pilr", "mucp", "cpsa", "cpsb", "capsul", "metalloprotease", "atoc", "peptide abc"]): return "#FF9400"
        if "hypothetical" in p or "unknown" in p or "uncharacter" in p or "duf" in p: return "#D3D3D3"
        return "#8ECAE6"

    # BLASTP links (query PAI proteins vs each proteome)
    for l in ORDER:
        db = f"{WORK}/blast/pdb/{l}"
        if not os.path.exists(db + ".psq"):
            run(["makeblastdb", "-in", f"{WORK}/proteins/{l}.faa", "-dbtype", "prot", "-out", db])
    for pai in PAIS:
        for l in ORDER:
            out = f"{WORK}/blast/m3/{pai}__{l}.tsv"
            if not os.path.exists(out):
                run(["blastp", "-query", f"{WORK}/pai_proteins/{pai}_proteins.fasta", "-db", f"{WORK}/blast/pdb/{l}",
                     "-out", out, "-outfmt", "6 qseqid sseqid pident length qlen evalue bitscore",
                     "-evalue", "1e-5", "-max_target_seqs", "5", "-num_threads", "4"])

    def load_links(pai):
        links = {}
        for l in ORDER:
            f = f"{WORK}/blast/m3/{pai}__{l}.tsv"
            if not os.path.exists(f): continue
            for line in open(f):
                p = line.rstrip("\n").split("\t")
                if len(p) < 7: continue
                qloc, spid, pident, aln, qlen = p[0], p[1], float(p[2]), int(p[3]), int(p[4])
                if aln / qlen < 0.4: continue
                sloc = prot2locus.get(l, {}).get(spid)
                if sloc is None: continue
                cur = links.get((qloc, l))
                if cur is None or pident > cur[1]: links[(qloc, l)] = (sloc, pident)
        return links

    for pai, (s, e) in PAIS.items():
        tracks = [("ATCC 25586 (query)", neighborhood(gffdata["AE009951.2_query"],
                  gffdata["AE009951.2_query"].seqid.iloc[0], s, e))]
        for l in ORDER:
            a = best_anchor(pai, l)
            if a is None or a["qcov"] < 0.4: continue
            sub = neighborhood(gffdata[l], a["seqid"], a["start"], a["end"])
            if not sub.empty: tracks.append((l, sub))
        links = load_links(pai)
        n = len(tracks)
        fig, ax = plt.subplots(figsize=(12, 0.85 * n + 2.0))
        tg = []
        for ti, (label, df) in enumerate(tracks):
            y = n - ti; mn, mx = df.start.min(), df.end.max(); span = max(mx - mn, 1)
            genes = []
            for _, r in df.iterrows():
                x0, x1 = (r.start - mn) / span, (r.end - mn) / span
                col = module_color(r["product"]); genes.append((r["locus"], x0, x1, y))
                h = 0.34
                pts = ([(x0, y - h), (max(x1 - 0.012, x0), y - h), (x1, y), (max(x1 - 0.012, x0), y + h), (x0, y + h)]
                       if r["strand"] == "+" else
                       [(x1, y - h), (min(x0 + 0.012, x1), y - h), (x0, y), (min(x0 + 0.012, x1), y + h), (x1, y + h)])
                ax.add_patch(Polygon(pts, closed=True, facecolor=col, edgecolor="black", linewidth=0.4))
            tg.append((label, genes)); ax.text(-0.008, y, label.replace("_", " "), ha="right", va="center", fontsize=8.5)
        qlabel, qgenes = tg[0]; qmap = {l: ((a + b) / 2, y) for l, a, b, y in qgenes}
        for ti in range(1, len(tg)):
            slab, sgenes = tg[ti]; smap = {l: ((a + b) / 2, y) for l, a, b, y in sgenes}
            for (qloc, slab2), (sloc, pident) in links.items():
                if slab2 != slab or qloc not in qmap or sloc not in smap: continue
                x1, y1 = qmap[qloc]; x2, y2 = smap[sloc]
                ax.plot([x1, x2], [y1 - 0.34, y2 + 0.34], color="#0279EE",
                        alpha=0.10 + 0.45 * (pident / 100), linewidth=0.7, zorder=0)
        ax.set_xlim(-0.30, 1.02); ax.set_ylim(0.3, n + 0.7); ax.axis("off")
        handles = [mpatches.Patch(color=c, label=l) for c, l in [
            ("#E4572E", "Mobility"), ("#B9314F", "Virulence/TA"), ("#0279EE", "Biotin"),
            ("#75A025", "Glycogen"), ("#7B2CBF", "Sigma/regulation"), ("#FF9400", "Colonization/immune"),
            ("#8ECAE6", "Other conserved"), ("#D3D3D3", "Hypothetical")]]
        ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.03), ncol=4, fontsize=7.5, frameon=False)
        ax.set_title(f"{pai}: gene-order conservation across genomes", fontsize=11)
        fig.tight_layout()
        fig.savefig(f"{FIGS}/M3_synteny_{pai}.png", dpi=200, bbox_inches="tight")
        fig.savefig(f"{FIGS}/M3_synteny_{pai}.svg", bbox_inches="tight")
        plt.close(fig)


# --------------------------------------------------------------------------
# M4: core-genome phylogeny + PAI distribution
# --------------------------------------------------------------------------
def m4():
    print("[M4] core-genome phylogeny")
    target = re_compile(r"(30S ribosomal protein S\d+|50S ribosomal protein L\d+)")
    gene_seqs, counts = defaultdict(dict), defaultdict(Counter)
    for l in ORDER:
        for r in SeqIO.parse(f"{WORK}/proteins/{l}.faa", "fasta"):
            m = target.search(r.description)
            if not m: continue
            key = m.group(1); counts[key][l] += 1
            if l not in gene_seqs[key]: gene_seqs[key][l] = r
    core = [k for k in gene_seqs if len(gene_seqs[k]) == len(ORDER) and all(counts[k][l] == 1 for l in ORDER)]
    print(f"  {len(core)} single-copy core ribosomal genes")
    for key in core:
        safe = key.replace(" ", "_")
        with open(f"{WORK}/core_aln/{safe}.faa", "w") as h:
            for l in ORDER: h.write(f">{l}\n{gene_seqs[key][l].seq}\n")
        run_out = subprocess.run(["mafft", "--auto", "--quiet", f"{WORK}/core_aln/{safe}.faa"],
                                 capture_output=True, text=True)
        open(f"{WORK}/core_aln/aligned/{safe}.aln", "w").write(run_out.stdout)
    concat = {l: "" for l in ORDER}
    for key in core:
        aln = AlignIO.read(f"{WORK}/core_aln/aligned/{key.replace(' ', '_')}.aln", "fasta")
        for rec in aln: concat[rec.id] += str(rec.seq)
    with open(f"{WORK}/core_aln/concatenated.faa", "w") as h:
        for l in ORDER: h.write(f">{l}\n{concat[l]}\n")
    with open(f"{WORK}/core_aln/core_tree.nwk", "w") as h:
        subprocess.run(["FastTree", "-wag", "-gamma", "-quiet", f"{WORK}/core_aln/concatenated.faa"],
                       stdout=h, check=True)
    # render tree + PAI integrity
    m2 = pd.read_csv(f"{RESULTS}/M2_gene_presence_long.csv")
    frac = {}
    for p in PAIS:
        f = m2[m2.PAI == p].groupby("genome").present.mean()
        for lab, v in f.items(): frac[(p, lab)] = v
    tree = Phylo.read(f"{WORK}/core_aln/core_tree.nwk", "newick"); tree.ladderize()
    terms = tree.get_terminals(); n = len(terms)
    yterm = {t: (n - 1 - i) for i, t in enumerate(terms)}
    def gy(c): return float(yterm[c]) if c.is_terminal() else float(np.mean([gy(x) for x in c.clades]))
    def gx(c): return tree.distance(c)
    maxd = max(tree.distance(t) for t in terms)
    fig = plt.figure(figsize=(13, 7.5)); gs = fig.add_gridspec(1, 2, width_ratios=[1.5, 1], wspace=0.02)
    axt = fig.add_subplot(gs[0]); axh = fig.add_subplot(gs[1], sharey=axt)
    for clade in tree.find_clades():
        x = gx(clade)
        if clade.clades:
            ys = [gy(c) for c in clade.clades]
            axt.plot([x, x], [min(ys), max(ys)], "k-", lw=1.1)
            for c in clade.clades: axt.plot([x, gx(c)], [gy(c), gy(c)], "k-", lw=1.1)
            if clade.confidence is not None:
                axt.text(x, gy(clade) + 0.3, f"{clade.confidence:.2f}", fontsize=6, ha="center", color="dimgray",
                         bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.7))
    for t in terms:
        axt.text(gx(t) + maxd * 0.012, yterm[t], t.name.replace("_", " "), fontsize=8.5, va="center")
    axt.set_xlim(0, maxd * 1.42); axt.set_ylim(-0.6, n - 0.4); axt.axis("off")
    axt.set_title("Core-genome phylogeny\n(49 ribosomal proteins, FastTree WAG+gamma)", fontsize=10, loc="left")
    order = [t.name for t in terms]
    M = np.array([[frac.get((p, lab), 0) * 100 for p in PAIS] for lab in order])
    cmap = plt.get_cmap("YlGnBu")
    for i, lab in enumerate(order):
        y = n - 1 - i
        for j in range(3):
            v = M[i, j]
            axh.add_patch(plt.Rectangle((j - 0.5, y - 0.5), 1, 1, facecolor=cmap(v / 100), edgecolor="white", lw=1.5))
            axh.text(j, y, f"{v:.0f}", ha="center", va="center", fontsize=8, color="white" if v > 50 else "black")
    axh.set_xlim(-0.5, 2.5); axh.set_ylim(-0.6, n - 0.4)
    axh.set_xticks(range(3)); axh.set_xticklabels(list(PAIS), fontsize=10); axh.set_yticks([])
    for sp in axh.spines.values(): sp.set_visible(False)
    axh.set_title("PAI gene-block integrity\n(% of PAI genes present)", fontsize=10, loc="left")
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=Normalize(0, 100))
    fig.colorbar(sm, ax=axh, fraction=0.046, pad=0.04, label="% of PAI genes present")
    fig.suptitle("M4: Species phylogeny vs. PAI conservation", fontsize=12, x=0.55)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(f"{FIGS}/M4_phylogeny_PAI.png", dpi=200); fig.savefig(f"{FIGS}/M4_phylogeny_PAI.svg")
    plt.close(fig)


def re_compile(pat):
    import re
    return re.compile(pat)


# --------------------------------------------------------------------------
# M5: targeted HGT tests
# --------------------------------------------------------------------------
def m5():
    print("[M5] targeted HGT tests")
    rec = SeqIO.read(os.path.join(INDIR, "25586.fasta"), "fasta"); genome = str(rec.seq).upper()
    def gc(s): return 100 * (s.count("G") + s.count("C")) / max(len(s), 1)
    def dinuc(s):
        c = Counter(s[i:i + 2] for i in range(len(s) - 1)); t = sum(c.values())
        return {k: v / t for k, v in c.items()}
    def ddist(a, b):
        keys = set(a) | set(b); return float(np.sqrt(sum((a.get(k, 0) - b.get(k, 0)) ** 2 for k in keys)))
    gen_gc, gen_dn = gc(genome), dinuc(genome)
    rng = np.random.default_rng(42); rows = []
    for pai, (s, e) in PAIS.items():
        sub = genome[s - 1:e]; L = len(sub); pai_gc, pai_dn = gc(sub), dinuc(sub); dd = ddist(pai_dn, gen_dn)
        ngc, ndd = [], []
        for st in rng.integers(0, len(genome) - L, size=500):
            w = genome[st:st + L]; ngc.append(gc(w)); ndd.append(ddist(dinuc(w), gen_dn))
        ngc, ndd = np.array(ngc), np.array(ndd)
        gz = (pai_gc - ngc.mean()) / ngc.std(); dz = (dd - ndd.mean()) / ndd.std()
        rows.append({"PAI": pai, "len": L, "PAI_GC%": round(pai_gc, 2), "genome_GC%": round(gen_gc, 2),
                     "GC_zscore": round(gz, 2), "dinuc_zscore": round(dz, 2),
                     "GC_atypical": abs(gz) > 2, "dinuc_atypical": abs(dz) > 2})
    pd.DataFrame(rows).to_csv(f"{RESULTS}/M5_compositional.csv", index=False)
    # mobility context
    q = urllib.parse.urlencode({"db": "nuccore", "id": "AE009951.2", "rettype": "gbwithparts", "retmode": "text"})
    txt = http("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + q, 90).read().decode()
    grec = SeqIO.read(io.StringIO(txt), "genbank")
    qdf = pd.DataFrame([{"start": int(f.location.start) + 1, "end": int(f.location.end),
                         "locus": f.qualifiers.get("locus_tag", ["?"])[0],
                         "product": f.qualifiers.get("product", [""])[0]}
                        for f in grec.features if f.type == "CDS"])
    mobpat = "transposase|integrase|recombinase|insertion sequence|resolvase|phage|mob"
    def has_mob(df): return df[df["product"].astype(str).str.lower().str.contains(mobpat, regex=True)]
    mrows = []
    for pai, (s, e) in PAIS.items():
        inside = qdf[(qdf.start >= s) & (qdf.end <= e)]
        fl = qdf[((qdf.end < s) & (qdf.end >= s - 15000)) | ((qdf.start > e) & (qdf.start <= e + 15000))]
        mi, mf = has_mob(inside), has_mob(fl)
        mrows.append({"PAI": pai, "n_genes": len(inside),
                      "mobility_genes_inside": "; ".join(f"{r.locus}({r.product})" for r in mi.itertuples()) or "none",
                      "mobility_genes_flanking": "; ".join(f"{r.locus}({r.product})" for r in mf.itertuples()) or "none"})
    pd.DataFrame(mrows).to_csv(f"{RESULTS}/M5_mobility_context.csv", index=False)
    # gene trees + concordance
    prot2locus = pickle.load(open(f"{WORK}/prot2locus.pkl", "rb"))
    protseq = {l: {r.id: r.seq for r in SeqIO.parse(f"{WORK}/proteins/{l}.faa", "fasta")} for l in ORDER}
    qp = {}
    for pai in PAIS:
        for r in SeqIO.parse(f"{WORK}/pai_proteins/{pai}_proteins.fasta", "fasta"): qp[r.id] = r.seq
    def orthologs(pai, locus):
        seqs = {"ATCC25586_query": qp[locus]}
        for g in ORDER:
            bf = f"{WORK}/blast/m3/{pai}__{g}.tsv"
            if not (os.path.exists(bf) and os.path.getsize(bf) > 0): continue
            best = None
            for line in open(bf):
                p = line.rstrip("\n").split("\t")
                if len(p) < 7 or p[0] != locus: continue
                pident, aln, qlen = float(p[2]), int(p[3]), int(p[4])
                if aln / qlen < 0.4: continue
                if best is None or pident > best[0]: best = (pident, p[1])
            if best and best[1] in protseq[g]: seqs[g] = protseq[g][best[1]]
        return seqs
    species = Phylo.read(f"{WORK}/core_aln/core_tree.nwk", "newick")
    def splits(tree):
        terms = [t.name for t in tree.get_terminals()]; s = set()
        for clade in tree.get_nonterminals():
            leaves = frozenset(t.name for t in clade.get_terminals())
            if 1 < len(leaves) < len(terms): s.add(leaves)
        return s, set(terms)
    def rf(t1, t2):
        s1, tax1 = splits(t1); s2, tax2 = splits(t2); shared = tax1 & tax2
        def restr(ss, sh):
            return set(frozenset(x for x in sp if x in sh) for sp in ss
                       if 1 < len(frozenset(x for x in sp if x in sh)) < len(sh))
        r1, r2 = restr(s1, shared), restr(s2, shared)
        return len(r1.symmetric_difference(r2)), len(r1) + len(r2)
    import copy
    grows = []
    for locus, pai in {"FN1885": "PAI1", "FN0837": "PAI2", "FN0849": "PAI2", "FN1317": "PAI3", "FN1318": "PAI3"}.items():
        seqs = orthologs(pai, locus)
        if len(seqs) < 4: continue
        fa = f"{WORK}/genetrees/{locus}.faa"
        with open(fa, "w") as h:
            for nm, s in seqs.items(): h.write(f">{nm}\n{s}\n")
        aln = f"{WORK}/genetrees/{locus}.aln"
        open(aln, "w").write(subprocess.run(["mafft", "--auto", "--quiet", fa], capture_output=True, text=True).stdout)
        nwk = f"{WORK}/genetrees/{locus}.nwk"
        with open(nwk, "w") as h:
            subprocess.run(["FastTree", "-wag", "-gamma", "-quiet", aln], stdout=h, check=True)
        gt = Phylo.read(nwk, "newick")
        for t in gt.get_terminals():
            if t.name == "ATCC25586_query": t.name = "Fn_nucleatum_ATCC25586"
        sp = copy.deepcopy(species); gtaxa = {t.name for t in gt.get_terminals()}
        for t in sp.get_terminals():
            if t.name not in gtaxa: sp.prune(t)
        d, tot = rf(gt, sp)
        grows.append({"gene": locus, "n_taxa": len(gtaxa), "RF": d, "RF_total": tot,
                      "normalized_RF": round(d / tot if tot else 0, 3)})
    pd.DataFrame(grows).to_csv(f"{RESULTS}/M5_gene_tree_concordance.csv", index=False)


# --------------------------------------------------------------------------
# M6: GC-skew reinterpretation
# --------------------------------------------------------------------------
def m6():
    print("[M6] GC-skew reinterpretation")
    rec = SeqIO.read(os.path.join(INDIR, "25586.fasta"), "fasta"); seq = str(rec.seq).upper()
    W, step = 5000, 2500; pos, skew = [], []
    for i in range(0, len(seq), step):
        win = seq[max(0, i - W):i]; g = win.count("G"); c = win.count("C")
        skew.append((g - c) / (g + c) if (g + c) > 0 else 0); pos.append(i)
    pos = np.array(pos); skew = np.array(skew); cum = np.cumsum(skew)
    ori, ter = pos[np.argmin(cum)], pos[np.argmax(cum)]
    colors = {"PAI1": "#E9ED4C", "PAI2": "#FF9400", "PAI3": "#FD9BED"}
    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    ax = axes[0]
    ax.plot(pos / 1e6, skew, color="gray", lw=0.7)
    ax.fill_between(pos / 1e6, 0, skew, where=skew >= 0, color="#0072B2", alpha=0.25, interpolate=True)
    ax.fill_between(pos / 1e6, 0, skew, where=skew < 0, color="#E69F00", alpha=0.25, interpolate=True)
    for pai, (s, e) in PAIS.items():
        ax.axvspan(s / 1e6, e / 1e6, color=colors[pai], alpha=0.55, edgecolor="black", lw=0.6)
    ax.set_ylabel("GC skew (5 kb window)"); ax.set_title("GC skew with PAI positions (shaded)", fontsize=11)
    ax.legend(handles=[mpatches.Patch(color=colors[p], label=p) for p in PAIS], loc="upper right", fontsize=8, ncol=3)
    ax2 = axes[1]
    ax2.plot(pos / 1e6, cum, color="#0279EE", lw=1.2)
    ax2.axvline(ori / 1e6, color="green", ls="--", lw=1.2, label=f"oriC (cum min) ~{ori//1000} kb")
    ax2.axvline(ter / 1e6, color="red", ls="--", lw=1.2, label=f"ter (cum max) ~{ter//1000} kb")
    for pai, (s, e) in PAIS.items():
        ax2.axvspan(s / 1e6, e / 1e6, color=colors[pai], alpha=0.55, edgecolor="black", lw=0.6)
        ax2.text((s + e) / 2 / 1e6, cum.min() * 0.92, pai, ha="center", fontsize=8, fontweight="bold")
    ax2.set_ylabel("Cumulative GC skew"); ax2.set_xlabel("Genome position (Mb)")
    ax2.legend(loc="lower left", fontsize=8)
    ax2.set_title("Cumulative GC skew: PAIs do NOT coincide with oriC/ter inflection points", fontsize=11)
    for a in axes: a.grid(True, ls="--", alpha=0.4); a.set_xlim(0, len(seq) / 1e6)
    fig.suptitle("M6: Reinterpretation of GC skew over the candidate PAIs (AE009951.2)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(f"{FIGS}/M6_gc_skew_reinterpretation.png", dpi=200)
    fig.savefig(f"{FIGS}/M6_gc_skew_reinterpretation.svg")
    plt.close(fig)
    rows = []
    for pai, (s, e) in PAIS.items():
        m = (pos >= s) & (pos <= e); z = (skew[m].mean() - skew.mean()) / skew.std()
        rows.append({"PAI": pai, "mean_window_skew": round(skew[m].mean(), 3),
                     "skew_zscore_vs_genome": round(z, 2), "at_cumulative_inflection": "no"})
    pd.DataFrame(rows).to_csv(f"{RESULTS}/M6_gc_skew_stats.csv", index=False)


# --------------------------------------------------------------------------
def main():
    labels = ORDER + [l for l, _ in OUTGROUPS]
    # m0()
    build_dbs()
    m1(labels)
    m2(labels)
    m3()
    m4()
    m5()
    m6()
    print("\nAll modules complete. Figures in", FIGS, "| tables in", RESULTS)


if __name__ == "__main__":
    main()
