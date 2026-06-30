#!/usr/bin/env python3

import json
import re
import shutil
import subprocess
import zipfile
from pathlib import Path
from datetime import datetime

# ------------------------------------------------------------
# Output folders
# ------------------------------------------------------------

OUTDIR = Path("ncbi_species_download")
FASTA_DIR = OUTDIR / "fastas"
ANNOT_DIR = OUTDIR / "annotations"
META_DIR = OUTDIR / "metadata"
ZIP_DIR = OUTDIR / "zips"

for d in [FASTA_DIR, ANNOT_DIR, META_DIR, ZIP_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------
# Species marked as N in your table
# ------------------------------------------------------------

SPECIES = [
    "Drosophila melanogaster",
    "Aedes aegypti",
    "Anopheles gambiae",
    "Ctenocephalides felis",
    "Bombyx mori",
    "Danaus plexippus",
    "Operophtera brumata",
    "Papilio xuthus",
    "Plutella xylostella",
    "Agrilus planipennis",
    "Nicrophorus vespilloides",
    "Onthophagus taurus",
    "Oryctes borbonicus",
    "Anoplophora glabripennis",
    "Leptotarsa decemlineata",
    "Diabrotica virgifera",
    "Dendroctonus ponderosae",
    "Aethina tumida",
    "Tribolium castaneum",
    "Asbolus verrucosus",
    "Apis mellifera",
    "Bombus impatiens",
    "Atta cephalotes",
    "Acromyrmex echinatior",
    "Harpegnathos saltator",
    "Solenopsis invicta",
    "Polistes dominula",
    "Polistes canadensis",
    "Nasonia vitripennis",
    "Cephus cinctus",
    "Orussus abietinus",
    "Athalia rosae",
    "Pediculus humanus",
    "Nilaparvata lugens",
    "Laodelphax striatellus",
    "Cimex lectularius",
    "Halyomorpha halys",
    "Bemisia tabaci",
    "Melanaphis sacchari",
    "Aphis gossypii",
    "Rhopalosiphum maidis",
    "Acyrthosiphon pisum",
    "Diuraphis noxia",
    "Myzus persicae",
    "Sipha flava",
    "Diaphorina citri",
    "Frankliniella occidentalis",
    "Zootermopsis nevadensis",
    "Cryptotermes secundus",
    "Blattella germanica",
    "Orchesella cincta",
    "Folsomia candida",
    "Daphnia pulex",
    "Tigriopus californicus",
    "Eurytemora affinis",
    "Hyalella azteca",
    "Armadillidium vulgare",
    "Penaeus vannamei",
    "Leptotrombidium deliense",
    "Dinothrombium tinctorium",
    "Tetranychus urticae",
    "Euroglyphus maynei",
    "Dermatophagoides pteronyssinus",
    "Sarcoptes scabiei",
    "Varroa destructor",
    "Varroa jacobsoni",
    "Tropilaelaps mercedesae",
    "Galendromus occidentalis",
    "Ixodes scapularis",
    "Trichonephila clavipes",
    "Parasteatoda tepidariorum",
    "Stegodyphus mimosarum",
    "Centruroides sculpturatus",
    "Limulus polyphemus",
    "Hypsibius dujardini",
    "Ramazzottius varieornatus",
    "Priapulus caudatus",
    "Nematostella vectensis",
    "Hydra vulgaris",
    "Trichoplax adhaerens",
    "Trichoplax H2",
    "Amphimedon queenslandica",
]

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

def safe_name(x):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", x.strip()).strip("_")


def run_cmd(cmd, allow_fail=False):
    print("[CMD]", " ".join(cmd))
    p = subprocess.run(cmd, text=True, capture_output=True)
    if p.returncode != 0 and not allow_fail:
        raise RuntimeError(
            f"Command failed:\n{' '.join(cmd)}\n\nSTDOUT:\n{p.stdout}\n\nSTDERR:\n{p.stderr}"
        )
    return p


def load_summary_json(stdout):
    """
    NCBI datasets normally returns JSON with a 'reports' field.
    This also tolerates JSONL-style output just in case.
    """
    stdout = stdout.strip()
    if not stdout:
        return []

    try:
        obj = json.loads(stdout)
        if isinstance(obj, dict) and "reports" in obj:
            return obj["reports"]
        if isinstance(obj, list):
            return obj
        return [obj]
    except json.JSONDecodeError:
        reports = []
        for line in stdout.splitlines():
            line = line.strip()
            if not line:
                continue
            try:
                reports.append(json.loads(line))
            except json.JSONDecodeError:
                pass
        return reports


def get_nested(d, paths, default=""):
    for path in paths:
        cur = d
        ok = True
        for key in path:
            if isinstance(cur, dict) and key in cur:
                cur = cur[key]
            else:
                ok = False
                break
        if ok and cur is not None:
            return cur
    return default


def has_annotation(report):
    annot = get_nested(
        report,
        [
            ["annotation_info"],
            ["annotationInfo"],
            ["assembly_info", "annotation_metadata"],
            ["assemblyInfo", "annotationMetadata"],
        ],
        default=None,
    )
    return annot not in [None, "", {}, []]


def score_report(report):
    """
    Pick the best assembly:
    1. RefSeq reference genome
    2. RefSeq representative genome
    3. Complete Genome
    4. Chromosome
    5. Scaffold
    6. Contig
    7. Has annotation
    8. Most recent release date
    """

    accession = get_nested(report, [["accession"]], "")
    refseq_category = str(get_nested(
        report,
        [["assembly_info", "refseq_category"], ["assemblyInfo", "refseqCategory"]],
        "",
    )).lower()

    assembly_level = str(get_nested(
        report,
        [["assembly_info", "assembly_level"], ["assemblyInfo", "assemblyLevel"]],
        "",
    )).lower()

    release_date = str(get_nested(
        report,
        [["assembly_info", "release_date"], ["assemblyInfo", "releaseDate"]],
        "",
    ))

    level_score = {
        "complete genome": 4,
        "chromosome": 3,
        "scaffold": 2,
        "contig": 1,
    }.get(assembly_level, 0)

    category_score = 0
    if "reference" in refseq_category:
        category_score = 3
    elif "representative" in refseq_category:
        category_score = 2

    source_score = 1 if accession.startswith("GCF_") else 0
    annot_score = 1 if has_annotation(report) else 0

    try:
        date_score = datetime.fromisoformat(release_date[:10]).timestamp()
    except Exception:
        date_score = 0

    return (
        category_score,
        source_score,
        level_score,
        annot_score,
        date_score,
    )


def classify_report(report):
    accession = get_nested(report, [["accession"]], "")

    organism = get_nested(
        report,
        [["organism", "organism_name"], ["organism", "organismName"]],
        "",
    )

    taxid = get_nested(
        report,
        [["organism", "tax_id"], ["organism", "taxId"]],
        "",
    )

    assembly_name = get_nested(
        report,
        [["assembly_info", "assembly_name"], ["assemblyInfo", "assemblyName"]],
        "",
    )

    assembly_level = get_nested(
        report,
        [["assembly_info", "assembly_level"], ["assemblyInfo", "assemblyLevel"]],
        "",
    )

    refseq_category = get_nested(
        report,
        [["assembly_info", "refseq_category"], ["assemblyInfo", "refseqCategory"]],
        "",
    )

    release_date = get_nested(
        report,
        [["assembly_info", "release_date"], ["assemblyInfo", "releaseDate"]],
        "",
    )

    submitter = get_nested(
        report,
        [["assembly_info", "submitter"], ["assemblyInfo", "submitter"]],
        "",
    )

    text = json.dumps(report).lower()

    # NCBI does not have a single universal "T2T = yes/no" field for all species.
    # This flags assemblies whose metadata explicitly says T2T / telomere-to-telomere / gapless.
    likely_t2t = bool(
        re.search(r"\bt2t\b|telomere[- ]to[- ]telomere|gapless", text, re.IGNORECASE)
    )

    level_lower = str(assembly_level).lower()

    return {
        "accession": accession,
        "organism_ncbi": organism,
        "taxid": taxid,
        "assembly_name": assembly_name,
        "assembly_level": assembly_level,
        "refseq_category": refseq_category,
        "release_date": release_date,
        "submitter": submitter,
        "is_complete_genome": level_lower == "complete genome",
        "is_chromosome_level": level_lower == "chromosome",
        "is_scaffold_level": level_lower == "scaffold",
        "is_contig_level": level_lower == "contig",
        "likely_t2t_from_metadata": likely_t2t,
        "annotation_available_in_metadata": has_annotation(report),
    }


def extract_relevant_files(zip_path, species, accession):
    sp_safe = safe_name(species)
    prefix = f"{sp_safe}__{accession}"

    fasta_count = 0
    annot_count = 0
    meta_count = 0

    with zipfile.ZipFile(zip_path, "r") as z:
        for member in z.namelist():
            lower = member.lower()
            filename = Path(member).name

            if not filename:
                continue

            # FASTA genome
            if lower.endswith(("_genomic.fna", "_genomic.fna.gz", ".fna", ".fna.gz")):
                out = FASTA_DIR / f"{prefix}__{filename}"
                with z.open(member) as src, open(out, "wb") as dst:
                    shutil.copyfileobj(src, dst)
                fasta_count += 1

            # Annotation
            elif lower.endswith((".gff", ".gff3", ".gtf", ".gbff", ".gff.gz", ".gff3.gz", ".gtf.gz", ".gbff.gz")):
                out = ANNOT_DIR / f"{prefix}__{filename}"
                with z.open(member) as src, open(out, "wb") as dst:
                    shutil.copyfileobj(src, dst)
                annot_count += 1

            # Metadata
            elif filename in ["assembly_data_report.jsonl", "dataset_catalog.json"]:
                out = META_DIR / f"{prefix}__{filename}"
                with z.open(member) as src, open(out, "wb") as dst:
                    shutil.copyfileobj(src, dst)
                meta_count += 1

    return fasta_count, annot_count, meta_count


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

rows = []

for species in SPECIES:
    print("\n" + "=" * 80)
    print(f"Searching NCBI for: {species}")
    print("=" * 80)

    sp_safe = safe_name(species)

    summary_cmd = [
        "datasets",
        "summary",
        "genome",
        "taxon",
        species,
    ]

    p = run_cmd(summary_cmd, allow_fail=True)

    if p.returncode != 0:
        rows.append({
            "query_species": species,
            "status": "summary_failed",
            "error": p.stderr.strip().replace("\n", " | "),
        })
        continue

    reports = load_summary_json(p.stdout)

    if not reports:
        rows.append({
            "query_species": species,
            "status": "no_assembly_found",
        })
        continue

    best = sorted(reports, key=score_report, reverse=True)[0]
    info = classify_report(best)
    accession = info["accession"]

    if not accession:
        rows.append({
            "query_species": species,
            "status": "no_accession_found",
        })
        continue

    print(f"Selected: {accession} | {info['assembly_level']} | {info['assembly_name']}")

    # Save raw selected metadata
    selected_meta = META_DIR / f"{sp_safe}__{accession}__selected_report.json"
    with open(selected_meta, "w") as f:
        json.dump(best, f, indent=2)

    zip_path = ZIP_DIR / f"{sp_safe}__{accession}.zip"

    download_cmd = [
        "datasets",
        "download",
        "genome",
        "accession",
        accession,
        "--include",
        "genome,gff3,gtf,gbff",
        "--filename",
        str(zip_path),
    ]

    p2 = run_cmd(download_cmd, allow_fail=True)

    if p2.returncode != 0:
        rows.append({
            "query_species": species,
            "status": "download_failed",
            "error": p2.stderr.strip().replace("\n", " | "),
            **info,
        })
        continue

    fasta_count, annot_count, meta_count = extract_relevant_files(zip_path, species, accession)

    rows.append({
        "query_species": species,
        "status": "downloaded",
        **info,
        "fasta_files_extracted": fasta_count,
        "annotation_files_extracted": annot_count,
        "metadata_files_extracted": meta_count,
        "zip_file": str(zip_path),
    })

# ------------------------------------------------------------
# Write reports
# ------------------------------------------------------------

report_tsv = OUTDIR / "ncbi_assembly_level_report.tsv"
report_json = OUTDIR / "ncbi_assembly_level_report.json"

all_keys = sorted(set().union(*(r.keys() for r in rows)))

with open(report_tsv, "w") as f:
    f.write("\t".join(all_keys) + "\n")
    for r in rows:
        f.write("\t".join(str(r.get(k, "")) for k in all_keys) + "\n")

with open(report_json, "w") as f:
    json.dump(rows, f, indent=2)

print("\nDone.")
print(f"FASTA folder:      {FASTA_DIR}")
print(f"Annotation folder: {ANNOT_DIR}")
print(f"Metadata folder:   {META_DIR}")
print(f"Report TSV:        {report_tsv}")
print(f"Report JSON:       {report_json}")

# Quick terminal summaries
def count_by(key):
    out = {}
    for r in rows:
        value = r.get(key, "NA")
        out[value] = out.get(value, 0) + 1
    return out

print("\nAssembly level summary:")
for k, v in count_by("assembly_level").items():
    print(f"  {k}: {v}")

print("\nLikely T2T from metadata:")
for r in rows:
    if str(r.get("likely_t2t_from_metadata", "")).lower() == "true":
        print(f"  {r.get('query_species')} | {r.get('accession')} | {r.get('assembly_name')}")
