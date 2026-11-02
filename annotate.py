#!/usr/bin/env python3
"""Clinical annotation of student_08.vcf (GRCh38).

Sources, queried on the date the script is run:
  * Ensembl VEP REST (consequences, MANE/canonical transcript choice)
  * ClinVar VCF GRCh38 (germline classification)
  * gnomAD GraphQL, dataset gnomad_r4 (joint exome+genome allele frequency)

The summary table keeps one row per input variant. The transcript is chosen
by Ensembl VEP --flag_pick with this order:
  mane_select, mane_plus_clinical, canonical, appris, tsl,
  biotype, ccds, rank, length.
The reported consequence is the most severe Sequence Ontology term on that
transcript. If VEP's variant-level most_severe_consequence differs, it is
kept in a separate column.
"""

from __future__ import annotations

import csv
import json
import subprocess
import time
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
VCF = ROOT / "student_08.vcf"
OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)

USER_AGENT = "nsu-bioinfo-hw4/student_08 (academic annotation)"
CLINVAR_VCF = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
GNOMAD_URL = "https://gnomad.broadinstitute.org/api"
GNOMAD_DATASET = "gnomad_r4"

VEP_PARAMS = {
    "canonical": "1",
    "mane": "1",
    "hgvs": "1",
    "numbers": "1",
    "protein": "1",
    "appris": "1",
    "ccds": "1",
    "uniprot": "1",
    "domains": "1",
    "xref_refseq": "1",
    "variant_class": "1",
    "flag_pick": "1",
    "pick_order": (
        "mane_select,mane_plus_clinical,canonical,appris,tsl,"
        "biotype,ccds,rank,length"
    ),
    "AlphaMissense": "1",
}

# Ensembl VEP severity rank. Lower index = more severe.
# https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html
SEVERITY = [
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "transcript_amplification",
    "feature_elongation",
    "feature_truncation",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",
    "splice_donor_5th_base_variant",
    "splice_region_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
    "incomplete_terminal_codon_variant",
    "start_retained_variant",
    "stop_retained_variant",
    "synonymous_variant",
    "coding_sequence_variant",
    "mature_miRNA_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
    "intron_variant",
    "NMD_transcript_variant",
    "non_coding_transcript_variant",
    "coding_transcript_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "TFBS_ablation",
    "TFBS_amplification",
    "TF_binding_site_variant",
    "regulatory_region_ablation",
    "regulatory_region_amplification",
    "regulatory_region_variant",
    "intergenic_variant",
]
SEVERITY_RANK = {term: i for i, term in enumerate(SEVERITY)}

TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}

PICK_RULE = (
    "VEP flag_pick: MANE Select > MANE Plus Clinical > canonical > "
    "APPRIS > TSL > biotype > CCDS > consequence rank > transcript length. "
    "On the chosen transcript the most severe SO term is reported."
)


def http_json(url: str, payload: dict | None = None, timeout: int = 180, retries: int = 4):
    data = None if payload is None else json.dumps(payload).encode()
    headers = {"User-Agent": USER_AGENT, "Accept": "application/json"}
    if payload is not None:
        headers["Content-Type"] = "application/json"
    last = None
    for attempt in range(retries):
        req = urllib.request.Request(url, data=data, headers=headers)
        try:
            with urllib.request.urlopen(req, timeout=timeout) as response:
                body = response.read().decode()
                return json.loads(body), dict(response.headers)
        except urllib.error.HTTPError as exc:
            detail = exc.read().decode(errors="replace")[:1500]
            last = f"HTTP {exc.code}: {detail}"
            if exc.code in (400, 404):
                break
        except Exception as exc:  # noqa: BLE001 - network retry
            last = f"{type(exc).__name__}: {exc}"
        time.sleep(2 + attempt * 3)
    raise RuntimeError(last or "request failed")


def read_vcf(path: Path) -> list[dict]:
    records = []
    with path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            chrom, pos, vid, ref, alt, qual, filt, info = line.rstrip("\n").split("\t")[:8]
            pair = (ref, alt)
            if pair in TRANSITIONS:
                kind = "transition"
            elif len(ref) == 1 and len(alt) == 1 and ref in "ACGT" and alt in "ACGT":
                kind = "transversion"
            else:
                kind = "other"
            records.append(
                {
                    "chrom": chrom,
                    "pos": int(pos),
                    "id": vid,
                    "ref": ref,
                    "alt": alt,
                    "qual": qual,
                    "filter": filt,
                    "info": info,
                    "substitution": kind,
                    "key": f"{chrom}:{pos}:{ref}:{alt}",
                    "vep_input": f"{chrom} {pos} {vid} {ref} {alt} . . .",
                    "gnomad_id": f"{chrom}-{pos}-{ref}-{alt}",
                }
            )
    return records


def write_tsv(path: Path, rows: list[dict], columns: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore", delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({col: row.get(col, "") for col in columns})


def most_severe(terms: list[str]) -> str:
    if not terms:
        return ""
    return min(terms, key=lambda term: SEVERITY_RANK.get(term, 10_000))


def fetch_ensembl_version() -> dict:
    url = "https://rest.ensembl.org/info/software?content-type=application/json"
    data, _ = http_json(url, timeout=60)
    return data


def fetch_vep(records: list[dict]) -> list[dict]:
    query = "&".join(f"{key}={value}" for key, value in VEP_PARAMS.items())
    url = f"https://rest.ensembl.org/vep/homo_sapiens/region?{query}"
    merged: list[dict] = []
    batch_size = 10
    for start in range(0, len(records), batch_size):
        batch = records[start : start + batch_size]
        payload = {"variants": [row["vep_input"] for row in batch]}
        print(f"VEP {start + 1}-{start + len(batch)} / {len(records)}", flush=True)
        try:
            data, _ = http_json(url, payload, timeout=300)
        except RuntimeError as exc:
            if "AlphaMissense" in str(exc) or "HTTP 400" in str(exc):
                print("VEP batch failed, retrying without AlphaMissense", flush=True)
                fallback = url.replace("&AlphaMissense=1", "")
                data, _ = http_json(fallback, payload, timeout=300)
            else:
                # Split a failed batch so one variant cannot drop the rest.
                print(f"VEP batch failed ({exc}); retrying one by one", flush=True)
                data = []
                for row in batch:
                    one = {"variants": [row["vep_input"]]}
                    try:
                        piece, _ = http_json(url, one, timeout=180)
                    except RuntimeError as one_exc:
                        print(f"  failed {row['id']}: {one_exc}", flush=True)
                        piece = []
                    data.extend(piece)
        merged.extend(data)
        time.sleep(1)
    return merged


GNOMAD_SELECTION = """
      variant_id
      chrom
      pos
      ref
      alt
      rsids
      exome { ac an ac_hom filters }
      genome { ac an ac_hom filters }
      joint {
        ac
        an
        hemizygote_count
        faf95 { popmax popmax_population }
      }
    """


def fetch_gnomad(records: list[dict]) -> dict:
    """Query gnomAD in small batches. The public API often returns 'Service overloaded'."""
    found: dict[str, dict] = {}
    batch_size = 5
    for start in range(0, len(records), batch_size):
        batch = records[start : start + batch_size]
        fields = []
        for index, row in enumerate(batch):
            fields.append(
                f'v{index}: variant(variantId: "{row["gnomad_id"]}", dataset: {GNOMAD_DATASET}) {{ {GNOMAD_SELECTION} }}'
            )
        query = "query {\n" + "\n".join(fields) + "\n}"
        print(f"gnomAD {start + 1}-{start + len(batch)} / {len(records)}", flush=True)
        data = None
        for attempt in range(6):
            try:
                data, _ = http_json(GNOMAD_URL, {"query": query}, timeout=120, retries=1)
            except RuntimeError as exc:
                print(f"  transport error attempt {attempt + 1}: {exc}", flush=True)
                time.sleep(5 + attempt * 5)
                continue
            if data.get("errors"):
                message = json.dumps(data["errors"])[:300]
                print(f"  API error attempt {attempt + 1}: {message}", flush=True)
                data = None
                time.sleep(8 + attempt * 7)
                continue
            break
        if not data or not data.get("data"):
            print("  batch failed, falling back to single variants", flush=True)
            for row in batch:
                single = (
                    "query { v0: variant(variantId: \""
                    + row["gnomad_id"]
                    + f"\", dataset: {GNOMAD_DATASET}) {{ {GNOMAD_SELECTION} }} }}"
                )
                one = None
                for attempt in range(5):
                    try:
                        one, _ = http_json(GNOMAD_URL, {"query": single}, timeout=90, retries=1)
                    except RuntimeError as exc:
                        print(f"  {row['id']} transport {exc}", flush=True)
                        time.sleep(4 + attempt * 4)
                        continue
                    if one.get("errors"):
                        print(f"  {row['id']} error {json.dumps(one['errors'])[:200]}", flush=True)
                        one = None
                        time.sleep(6 + attempt * 5)
                        continue
                    break
                found[row["key"]] = (one or {}).get("data", {}).get("v0") if one else None
                time.sleep(1.5)
            continue
        payload = data["data"]
        for index, row in enumerate(batch):
            found[row["key"]] = payload.get(f"v{index}")
        time.sleep(2)
    return found


def fetch_gnomad_populations(variant_ids: list[str]) -> dict[str, list]:
    """Ancestry counts for a few variants discussed in the report."""
    out = {}
    selection = "joint { populations { id ac an } }"
    for variant_id in variant_ids:
        query = (
            "query { v0: variant(variantId: \""
            + variant_id
            + f"\", dataset: {GNOMAD_DATASET}) {{ {selection} }} }}"
        )
        data, _ = http_json(GNOMAD_URL, {"query": query}, timeout=90)
        joint = ((data.get("data") or {}).get("v0") or {}).get("joint") or {}
        out[variant_id] = joint.get("populations") or []
        time.sleep(1)
    return out


def fetch_clinvar(records: list[dict]) -> tuple[str, list[str]]:
    regions = [f"{row['chrom']}:{row['pos']}-{row['pos']}" for row in records]
    print(f"ClinVar tabix, {len(regions)} sites", flush=True)
    header = subprocess.run(
        ["tabix", "-H", CLINVAR_VCF],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    body = subprocess.run(
        ["tabix", CLINVAR_VCF, *regions],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    return header, body.splitlines()


def parse_info(info: str) -> dict[str, str]:
    parsed = {}
    for item in info.split(";"):
        if "=" not in item:
            parsed[item] = "1"
            continue
        key, value = item.split("=", 1)
        parsed[key] = value
    return parsed


def af_text(ac, an) -> str:
    if ac is None or an in (None, 0):
        return ""
    return f"{ac / an:.6g}"


def choose_transcript(consequences: list[dict]) -> dict | None:
    """Use the VEP PICK flag. Fall back to the same priority if the flag is absent."""
    if not consequences:
        return None
    picked = [row for row in consequences if row.get("pick") == 1]
    if len(picked) == 1:
        return picked[0]
    if len(picked) > 1:
        pool = picked
    else:
        pool = consequences

    def sort_key(row: dict):
        terms = row.get("consequence_terms") or []
        return (
            0 if row.get("mane_select") else 1,
            0 if row.get("mane_plus_clinical") else 1,
            0 if row.get("canonical") == 1 else 1,
            0 if str(row.get("biotype", "")).startswith("protein_coding") else 1,
            SEVERITY_RANK.get(most_severe(terms), 10_000),
            row.get("transcript_id") or "",
        )

    return sorted(pool, key=sort_key)[0]


def domain_text(row: dict) -> str:
    domains = row.get("domains") or []
    names = []
    for domain in domains:
        name = domain.get("name") or domain.get("description") or ""
        db = domain.get("db") or ""
        if name:
            names.append(f"{db}:{name}" if db else name)
    return "; ".join(names)


def consequence_row(vep: dict) -> dict:
    transcripts = vep.get("transcript_consequences") or []
    chosen = choose_transcript(transcripts)
    genes = sorted({row.get("gene_symbol") or row.get("gene_id") or "" for row in transcripts if row.get("gene_symbol") or row.get("gene_id")})
    if chosen is None:
        return {
            "gene": "",
            "all_genes": "; ".join(genes),
            "n_transcripts": 0,
            "transcript_id": "",
            "mane_select": "",
            "canonical": "",
            "biotype": "",
            "consequence": vep.get("most_severe_consequence") or "",
            "all_consequences_on_transcript": "",
            "impact": "",
            "hgvsc": "",
            "hgvsp": "",
            "amino_acids": "",
            "codons": "",
            "exon": "",
            "intron": "",
            "protein_position": "",
            "sift": "",
            "polyphen": "",
            "alphamissense": "",
            "appris": "",
            "ccds": "",
            "refseq": "",
            "domains": "",
            "pick_source": "none",
        }
    terms = chosen.get("consequence_terms") or []
    am_class = chosen.get("am_class") or chosen.get("alphamissense_class") or ""
    am_score = chosen.get("am_pathogenicity") or chosen.get("alphamissense_pathogenicity")
    if am_score is None and isinstance(chosen.get("alphamissense"), dict):
        am_score = chosen["alphamissense"].get("am_pathogenicity")
        am_class = am_class or chosen["alphamissense"].get("am_class", "")
    am = ""
    if am_score is not None or am_class:
        am = f"{am_class} {am_score}".strip()
    sift = ""
    if chosen.get("sift_prediction"):
        sift = f"{chosen.get('sift_prediction')} ({chosen.get('sift_score')})"
    poly = ""
    if chosen.get("polyphen_prediction"):
        poly = f"{chosen.get('polyphen_prediction')} ({chosen.get('polyphen_score')})"
    return {
        "gene": chosen.get("gene_symbol") or chosen.get("gene_id") or "",
        "all_genes": "; ".join(genes),
        "n_transcripts": len(transcripts),
        "transcript_id": chosen.get("transcript_id") or "",
        "mane_select": chosen.get("mane_select") or "",
        "canonical": "yes" if chosen.get("canonical") == 1 else "no",
        "biotype": chosen.get("biotype") or "",
        "consequence": most_severe(terms),
        "all_consequences_on_transcript": "; ".join(terms),
        "impact": chosen.get("impact") or "",
        "hgvsc": chosen.get("hgvsc") or "",
        "hgvsp": chosen.get("hgvsp") or "",
        "amino_acids": chosen.get("amino_acids") or "",
        "codons": chosen.get("codons") or "",
        "exon": chosen.get("exon") or "",
        "intron": chosen.get("intron") or "",
        "protein_position": chosen.get("protein_start") or "",
        "sift": sift,
        "polyphen": poly,
        "alphamissense": am,
        "appris": chosen.get("appris") or "",
        "ccds": chosen.get("ccds") or "",
        "refseq": "; ".join(chosen.get("refseq_transcript_ids") or []) or chosen.get("mane_select") or "",
        "domains": domain_text(chosen),
        "pick_source": "vep_pick" if chosen.get("pick") == 1 else "fallback_same_priority",
    }


def index_vep(vep_rows: list[dict], records: list[dict]) -> dict[str, dict]:
    by_input = {row.get("input"): row for row in vep_rows}
    by_id = {}
    for row in vep_rows:
        if row.get("id"):
            by_id.setdefault(row["id"], row)
    indexed = {}
    for record in records:
        hit = by_input.get(record["vep_input"]) or by_id.get(record["id"])
        if hit is None:
            # Coordinate fallback: start + allele string.
            allele = f"{record['ref']}/{record['alt']}"
            for row in vep_rows:
                if (
                    str(row.get("seq_region_name")) == record["chrom"]
                    and int(row.get("start", -1)) == record["pos"]
                    and row.get("allele_string") == allele
                ):
                    hit = row
                    break
        indexed[record["key"]] = hit or {}
    return indexed


def index_clinvar(lines: list[str], records: list[dict]) -> dict[str, list[dict]]:
    wanted = {(row["chrom"], row["pos"], row["ref"], row["alt"]): row["key"] for row in records}
    grouped: dict[str, list[dict]] = {row["key"]: [] for row in records}
    for line in lines:
        if not line or line.startswith("#"):
            continue
        chrom, pos, vid, ref, alt, qual, filt, info = line.split("\t")[:8]
        alts = alt.split(",")
        for allele in alts:
            key = wanted.get((chrom, int(pos), ref, allele))
            if key is None:
                continue
            parsed = parse_info(info)
            grouped[key].append(
                {
                    "clinvar_variation_id": vid,
                    "ref": ref,
                    "alt": allele,
                    "allele_id": parsed.get("ALLELEID", ""),
                    "clinsig": parsed.get("CLNSIG", "").replace("_", " "),
                    "clinsig_raw": parsed.get("CLNSIG", ""),
                    "clinsig_conflict": parsed.get("CLNSIGCONF", "").replace("_", " "),
                    "review_status": parsed.get("CLNREVSTAT", "").replace("_", " "),
                    "condition": parsed.get("CLNDN", "").replace("_", " ").replace("|", "; "),
                    "condition_ids": parsed.get("CLNDISDB", ""),
                    "hgvs": parsed.get("CLNHGVS", ""),
                    "molecular_consequence": parsed.get("MC", "").replace("|", " "),
                    "geneinfo": parsed.get("GENEINFO", ""),
                    "rs": parsed.get("RS", ""),
                    "origin": parsed.get("ORIGIN", ""),
                    "oncogenicity": parsed.get("ONC", ""),
                    "clinical_impact": parsed.get("SCISIG", ""),
                }
            )
    return grouped


def flatten_vep(records: list[dict], vep_by_key: dict[str, dict]) -> list[dict]:
    rows = []
    for record in records:
        vep = vep_by_key.get(record["key"]) or {}
        for transcript in vep.get("transcript_consequences") or []:
            terms = transcript.get("consequence_terms") or []
            rows.append(
                {
                    "chrom": record["chrom"],
                    "pos": record["pos"],
                    "rsid": record["id"],
                    "ref": record["ref"],
                    "alt": record["alt"],
                    "variant_class": vep.get("variant_class") or "",
                    "most_severe_consequence": vep.get("most_severe_consequence") or "",
                    "gene_symbol": transcript.get("gene_symbol") or "",
                    "gene_id": transcript.get("gene_id") or "",
                    "transcript_id": transcript.get("transcript_id") or "",
                    "biotype": transcript.get("biotype") or "",
                    "consequence_terms": "; ".join(terms),
                    "most_severe_on_transcript": most_severe(terms),
                    "impact": transcript.get("impact") or "",
                    "canonical": transcript.get("canonical") or "",
                    "mane_select": transcript.get("mane_select") or "",
                    "mane_plus_clinical": transcript.get("mane_plus_clinical") or "",
                    "pick": transcript.get("pick") or "",
                    "appris": transcript.get("appris") or "",
                    "ccds": transcript.get("ccds") or "",
                    "hgvsc": transcript.get("hgvsc") or "",
                    "hgvsp": transcript.get("hgvsp") or "",
                    "amino_acids": transcript.get("amino_acids") or "",
                    "codons": transcript.get("codons") or "",
                    "exon": transcript.get("exon") or "",
                    "intron": transcript.get("intron") or "",
                    "sift_prediction": transcript.get("sift_prediction") or "",
                    "sift_score": transcript.get("sift_score") if transcript.get("sift_score") is not None else "",
                    "polyphen_prediction": transcript.get("polyphen_prediction") or "",
                    "polyphen_score": transcript.get("polyphen_score") if transcript.get("polyphen_score") is not None else "",
                    "am_class": transcript.get("am_class") or "",
                    "am_pathogenicity": transcript.get("am_pathogenicity") if transcript.get("am_pathogenicity") is not None else "",
                    "strand": transcript.get("strand") or "",
                }
            )
    return rows


def build_summary(records, vep_by_key, gnomad, clinvar) -> list[dict]:
    rows = []
    for record in records:
        vep = vep_by_key.get(record["key"]) or {}
        chosen = consequence_row(vep)
        g = gnomad.get(record["key"]) or {}
        exome = (g or {}).get("exome") or {}
        genome = (g or {}).get("genome") or {}
        joint = (g or {}).get("joint") or {}
        exome_filters = exome.get("filters") or []
        genome_filters = genome.get("filters") or []
        if not g:
            joint_af = ""
            gnomad_note = "нет записи"
        elif joint.get("ac") in (0, None) and ("AC0" in exome_filters or "AC0" in genome_filters):
            # AC0 means the allele was not retained in the filtered gnomAD call set.
            joint_af = "нет в отфильтрованном наборе (AC0)"
            gnomad_note = "запись есть, AC=0 после фильтра AC0"
        else:
            joint_af = af_text(joint.get("ac"), joint.get("an"))
            notes = []
            if exome_filters:
                notes.append("exome:" + ",".join(exome_filters))
            if genome_filters:
                notes.append("genome:" + ",".join(genome_filters))
            gnomad_note = "; ".join(notes)
        faf = (joint.get("faf95") or {}) if joint else {}
        hits = clinvar.get(record["key"]) or []
        if not hits:
            present = "нет"
            sig = ""
            review = ""
            condition = ""
            conflict = ""
            vid = ""
            hgvs = ""
            mc = ""
        else:
            present = "да"
            sigs = sorted({hit["clinsig"] for hit in hits if hit["clinsig"]})
            sig = " | ".join(sigs) if sigs else "агрегатная классификация отсутствует"
            review = " | ".join(sorted({hit["review_status"] for hit in hits if hit["review_status"]}))
            condition = " | ".join(sorted({hit["condition"] for hit in hits if hit["condition"]}))
            conflict = " | ".join(sorted({hit["clinsig_conflict"] for hit in hits if hit["clinsig_conflict"]}))
            vid = "; ".join(hit["clinvar_variation_id"] for hit in hits)
            hgvs = "; ".join(hit["hgvs"] for hit in hits if hit["hgvs"])
            mc = "; ".join(hit["molecular_consequence"] for hit in hits if hit["molecular_consequence"])
        most = vep.get("most_severe_consequence") or ""
        differs = "yes" if most and chosen["consequence"] and most != chosen["consequence"] else "no"
        rows.append(
            {
                "rsid": record["id"],
                "coordinates": f"{record['chrom']}:{record['pos']} {record['ref']}>{record['alt']}",
                "chrom": record["chrom"],
                "pos": record["pos"],
                "ref": record["ref"],
                "alt": record["alt"],
                "substitution": record["substitution"],
                "gene": chosen["gene"],
                "all_genes": chosen["all_genes"],
                "n_transcripts": chosen["n_transcripts"],
                "transcript_id": chosen["transcript_id"],
                "transcript_rule": PICK_RULE,
                "mane_select": chosen["mane_select"],
                "canonical": chosen["canonical"],
                "biotype": chosen["biotype"],
                "consequence": chosen["consequence"],
                "impact": chosen["impact"],
                "all_consequences_on_transcript": chosen["all_consequences_on_transcript"],
                "vep_most_severe_consequence": most,
                "most_severe_differs_from_pick": differs,
                "hgvsc": chosen["hgvsc"],
                "hgvsp": chosen["hgvsp"],
                "amino_acids": chosen["amino_acids"],
                "codons": chosen["codons"],
                "exon": chosen["exon"],
                "intron": chosen["intron"],
                "protein_position": chosen["protein_position"],
                "sift": chosen["sift"],
                "polyphen": chosen["polyphen"],
                "alphamissense": chosen["alphamissense"],
                "appris": chosen["appris"],
                "ccds": chosen["ccds"],
                "domains": chosen["domains"],
                "variant_class": vep.get("variant_class") or "",
                "gnomad_af": joint_af,
                "gnomad_qc_note": gnomad_note,
                "gnomad_ac": joint.get("ac") if g else "",
                "gnomad_an": joint.get("an") if g else "",
                "gnomad_exome_af": af_text(exome.get("ac"), exome.get("an")) if exome else "",
                "gnomad_exome_ac": exome.get("ac") if exome else "",
                "gnomad_exome_an": exome.get("an") if exome else "",
                "gnomad_genome_af": af_text(genome.get("ac"), genome.get("an")) if genome else "",
                "gnomad_genome_ac": genome.get("ac") if genome else "",
                "gnomad_genome_an": genome.get("an") if genome else "",
                "gnomad_faf95_popmax": faf.get("popmax") if faf else "",
                "gnomad_faf95_popmax_population": faf.get("popmax_population") if faf else "",
                "gnomad_exome_filters": "; ".join(exome.get("filters") or []) if exome else "",
                "gnomad_genome_filters": "; ".join(genome.get("filters") or []) if genome else "",
                "gnomad_present": "да" if g else "нет",
                "clinvar_record": present,
                "clinvar_significance": sig,
                "clinvar_conflict_detail": conflict,
                "clinvar_review_status": review,
                "clinvar_condition": condition,
                "clinvar_variation_id": vid,
                "clinvar_hgvs": hgvs,
                "clinvar_molecular_consequence": mc,
                "pick_source": chosen["pick_source"],
            }
        )
    return rows


def gnomad_rows(records, gnomad) -> list[dict]:
    rows = []
    for record in records:
        g = gnomad.get(record["key"])
        if not g:
            rows.append(
                {
                    "rsid": record["id"],
                    "gnomad_id": record["gnomad_id"],
                    "present": "нет",
                    "dataset": GNOMAD_DATASET,
                }
            )
            continue
        exome = g.get("exome") or {}
        genome = g.get("genome") or {}
        joint = g.get("joint") or {}
        faf = joint.get("faf95") or {}
        pops = joint.get("populations") or []
        pop_af = []
        for pop in pops:
            ac, an = pop.get("ac"), pop.get("an")
            if an:
                pop_af.append(f"{pop.get('id')}={ac}/{an} ({ac / an:.6g})")
        rows.append(
            {
                "rsid": record["id"],
                "gnomad_id": g.get("variant_id") or record["gnomad_id"],
                "present": "да",
                "dataset": GNOMAD_DATASET,
                "rsids_gnomad": "; ".join(g.get("rsids") or []),
                "joint_ac": joint.get("ac"),
                "joint_an": joint.get("an"),
                "joint_af": af_text(joint.get("ac"), joint.get("an")),
                "joint_hemizygote_count": joint.get("hemizygote_count"),
                "exome_ac": exome.get("ac", ""),
                "exome_an": exome.get("an", ""),
                "exome_af": af_text(exome.get("ac"), exome.get("an")) if exome else "",
                "exome_ac_hom": exome.get("ac_hom", ""),
                "exome_filters": "; ".join(exome.get("filters") or []),
                "genome_ac": genome.get("ac", ""),
                "genome_an": genome.get("an", ""),
                "genome_af": af_text(genome.get("ac"), genome.get("an")) if genome else "",
                "genome_ac_hom": genome.get("ac_hom", ""),
                "genome_filters": "; ".join(genome.get("filters") or []),
                "faf95_popmax": faf.get("popmax", ""),
                "faf95_popmax_population": faf.get("popmax_population", ""),
                "joint_populations": "; ".join(pop_af),
            }
        )
    return rows


def clinvar_rows(records, clinvar) -> list[dict]:
    rows = []
    for record in records:
        hits = clinvar.get(record["key"]) or []
        if not hits:
            rows.append(
                {
                    "rsid": record["id"],
                    "chrom": record["chrom"],
                    "pos": record["pos"],
                    "ref": record["ref"],
                    "alt": record["alt"],
                    "clinvar_record": "нет",
                }
            )
            continue
        for hit in hits:
            rows.append(
                {
                    "rsid": record["id"],
                    "chrom": record["chrom"],
                    "pos": record["pos"],
                    "ref": record["ref"],
                    "alt": record["alt"],
                    "clinvar_record": "да",
                    **hit,
                }
            )
    return rows


def main() -> None:
    records = read_vcf(VCF)
    n21 = sum(row["chrom"] == "21" for row in records)
    n22 = sum(row["chrom"] == "22" for row in records)
    n_ti = sum(row["substitution"] == "transition" for row in records)
    n_tv = sum(row["substitution"] == "transversion" for row in records)
    stats = {
        "file": VCF.name,
        "n_records": len(records),
        "n_snv": sum(len(row["ref"]) == 1 and len(row["alt"]) == 1 for row in records),
        "chrom_21": n21,
        "chrom_22": n22,
        "transitions": n_ti,
        "transversions": n_tv,
        "ti_tv": round(n_ti / n_tv, 3) if n_tv else None,
        "assembly": "GRCh38",
    }
    (OUT / "vcf_statistics.json").write_text(json.dumps(stats, ensure_ascii=False, indent=2))
    print(json.dumps(stats, ensure_ascii=False), flush=True)

    version = fetch_ensembl_version()
    (OUT / "ensembl_software.json").write_text(json.dumps(version, indent=2))

    vep_path = OUT / "vep_raw.json"
    if vep_path.exists() and vep_path.stat().st_size > 1000:
        print("reusing existing VEP JSON", flush=True)
        vep = json.loads(vep_path.read_text())
    else:
        vep = fetch_vep(records)
        vep_path.write_text(json.dumps(vep, ensure_ascii=False))
    print(f"VEP records {len(vep)}", flush=True)

    gnomad_path = OUT / "gnomad_raw.json"
    if gnomad_path.exists() and gnomad_path.stat().st_size > 1000:
        print("reusing existing gnomAD JSON", flush=True)
        gnomad = json.loads(gnomad_path.read_text())
    else:
        gnomad = fetch_gnomad(records)
        gnomad_path.write_text(json.dumps(gnomad, ensure_ascii=False))

    header, lines = fetch_clinvar(records)
    (OUT / "clinvar_matches.vcf").write_text(header + "".join(line + "\n" for line in lines))
    print(f"ClinVar lines {len(lines)}", flush=True)

    vep_by_key = index_vep(vep, records)
    missing_vep = [row["id"] for row in records if not vep_by_key.get(row["key"])]
    if missing_vep:
        print("WARNING missing VEP", missing_vep, flush=True)
    clinvar = index_clinvar(lines, records)

    transcript_rows = flatten_vep(records, vep_by_key)
    transcript_cols = list(transcript_rows[0].keys()) if transcript_rows else ["rsid"]
    write_tsv(OUT / "vep_transcripts.tsv", transcript_rows, transcript_cols)

    g_rows = gnomad_rows(records, gnomad)
    write_tsv(OUT / "gnomad_annotation.tsv", g_rows, list(g_rows[0].keys()))

    c_rows = clinvar_rows(records, clinvar)
    c_cols = []
    for row in c_rows:
        for key in row:
            if key not in c_cols:
                c_cols.append(key)
    write_tsv(OUT / "clinvar_annotation.tsv", c_rows, c_cols)

    summary = build_summary(records, vep_by_key, gnomad, clinvar)
    columns = list(summary[0].keys())
    write_tsv(OUT / "annotated_variants.tsv", summary, columns)
    with (OUT / "annotated_variants.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(summary)

    russian_columns = [
        "идентификатор",
        "координаты_GRCh38",
        "ген",
        "другие_гены_VEP",
        "последствие_VEP",
        "самое_тяжёлое_последствие_по_всем_транскриптам",
        "транскрипт",
        "MANE_Select",
        "правило_выбора_транскрипта",
        "rsID",
        "gnomAD_AF",
        "запись_ClinVar",
        "интерпретация_ClinVar",
        "статус_ревью_ClinVar",
        "фенотип_ClinVar",
    ]
    russian_rows = [
        {
            "идентификатор": row["rsid"],
            "координаты_GRCh38": row["coordinates"],
            "ген": row["gene"],
            "другие_гены_VEP": row["all_genes"],
            "последствие_VEP": row["consequence"],
            "самое_тяжёлое_последствие_по_всем_транскриптам": row["vep_most_severe_consequence"],
            "транскрипт": row["transcript_id"],
            "MANE_Select": row["mane_select"],
            "правило_выбора_транскрипта": row["transcript_rule"],
            "rsID": row["rsid"],
            "gnomAD_AF": row["gnomad_af"],
            "запись_ClinVar": row["clinvar_record"],
            "интерпретация_ClinVar": row["clinvar_significance"],
            "статус_ревью_ClinVar": row["clinvar_review_status"],
            "фенотип_ClinVar": row["clinvar_condition"],
        }
        for row in summary
    ]
    write_tsv(OUT / "annotated_variants_ru.tsv", russian_rows, russian_columns)

    settings = {
        "input_vcf": VCF.name,
        "assembly": "GRCh38",
        "n_variants": len(records),
        "queried_on": time.strftime("%Y-%m-%d"),
        "ensembl_software": version,
        "vep_endpoint": "POST https://rest.ensembl.org/vep/homo_sapiens/region",
        "vep_parameters": VEP_PARAMS,
        "transcript_rule": PICK_RULE,
        "clinvar": {
            "source": CLINVAR_VCF,
            "access": "tabix query of each variant position, allele matched on REF and ALT",
            "file_date_header": "see clinvar_matches.vcf ##fileDate",
            "classification_field": "CLNSIG (aggregate germline classification); conflicts in CLNSIGCONF",
        },
        "gnomad": {
            "endpoint": GNOMAD_URL,
            "dataset": GNOMAD_DATASET,
            "assembly": "GRCh38",
            "allele_frequency": "joint AC/AN (exomes + genomes). Exome-only and genome-only AF are stored separately.",
            "faf95": "filtering allele frequency, 95% confidence, population maximum",
        },
        "summary_table": "results/annotated_variants.tsv (one row per input variant)",
        "vcf_statistics": stats,
    }
    (OUT / "annotation_settings.json").write_text(json.dumps(settings, ensure_ascii=False, indent=2))
    print(f"summary rows {len(summary)}", flush=True)
    print(f"transcript rows {len(transcript_rows)}", flush=True)
    print(f"clinvar present {sum(row['clinvar_record']=='да' for row in summary)}", flush=True)
    print(f"gnomad present {sum(row['gnomad_present']=='да' for row in summary)}", flush=True)


if __name__ == "__main__":
    main()
