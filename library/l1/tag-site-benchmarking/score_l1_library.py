#!/usr/bin/env python3
"""Score every L1 library proteoform on its own sequence, and publish to D1.

L1 is the 7,066-candidate library in `candidate_manifest.tsv` — every
`candidate_id` reads `l1-pNNNN`. Its `key` names the proteoform three ways:

    lib:ABCA8                 library_reference, the human canonical form
    iso:ALCAM:Q13740-2        an alternative isoform, by its own accession
    ort:cyno:MMP16:G7PC63     an ortholog, by its own accession

Everything published so far scored the **canonical** form. For the 47 library
members whose construct is not canonical that describes a different protein, and
for the 2,817 ortholog constructs it describes a different species. This scores
each proteoform from its own stored topology and disorder, keyed on
`uniprot_acc_full`, so the residue numbers come back in the construct's own
coordinates.

**Three of the features are canonical-coordinate.** ProDomino, EpicTope and
INDELi-E ran on the literature cohort in canonical numbering, so they are used
only where the proteoform IS the canonical form; everywhere else they are left
missing, which the logistic imputes. They are worth roughly 0.07 AUROC on the
target-class metric, so an isoform or ortholog score is the weaker claim and the
output says which rows carried them.

**The sequence is verified, not the length.** Two proteoforms of equal length
can still be different proteins — an allele variant, a different transcript of
the same size — so each L1 ORF is translated from `candidate_cds.json` and
compared residue by residue against the stored proteoform's sequence. Anything
that is not an exact match is reported with its first mismatch and skipped,
rather than scored against a sequence it is not.

    uv run --locked python score_l1_library.py --manifest <candidate_manifest.tsv>
    uv run --locked python score_l1_library.py --manifest ... --write-d1
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
import urllib.request
import warnings
from collections import defaultdict
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.Seq import Seq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[3] / "src"))

from train_site_models import FUSED_NUMERIC, models  # noqa: E402

PREDICTOR = "deliverome_fused_logistic"
PREDICTOR_VERSION = "1.0"
# Canonical-coordinate features: valid only when the proteoform is canonical.
CANONICAL_ONLY = [
    "prodomino",
    "prodomino_w5",
    "epictope",
    "epictope_w5",
    "indeli_e",
    "indeli_e_w5",
    "lane_disorder",
    "lane_surface_loop",
    "lane_terminal",
    "lane_snorkel",
    "lane_conservation",
    "contact_here",
    "dist_to_contact",
    "n_contact_w10",
    "contact_ecd_here",
    "dist_to_contact_ecd",
    "dist_to_pocket",
    "cons_entropy",
    "cons_entropy_w5",
    "cons_entropy_w15",
    "msa_depth",
]


def d1(sql: str, params=None, tries: int = 4) -> list[dict]:
    # saved D1 copy by default; live D1 only with S5E13_D1_LIVE=1 (d1_snapshot.py)
    import d1_snapshot

    if not d1_snapshot.live():
        return d1_snapshot.query(sql, params)
    acct = os.environ["CLOUDFLARE_ACCOUNT_ID"]
    token = os.environ["CLOUDFLARE_API_TOKEN"]
    db = os.environ["CLOUDFLARE_D1_SURFACEOME_PUBLIC_ID"]
    body = json.dumps({"sql": sql, "params": params or []}).encode()
    for attempt in range(tries):
        try:
            req = urllib.request.Request(
                f"https://api.cloudflare.com/client/v4/accounts/{acct}/d1/database/{db}/query",
                data=body,
                headers={"Authorization": f"Bearer {token}", "Content-Type": "application/json"},
            )
            with urllib.request.urlopen(req, timeout=180) as r:
                p = json.loads(r.read().decode())
            if not p.get("success"):
                raise RuntimeError(p.get("errors"))
            return p["result"][0]["results"]
        except Exception:
            if attempt == tries - 1:
                raise
            time.sleep(2 * (attempt + 1))
    return []


def unit(h):
    return [b / 255.0 for b in bytes.fromhex(h)] if h else None


def plddt(h):
    return [b / 2.55 for b in bytes.fromhex(h)] if h else None


def at(arr, res):
    if arr is None:
        return None
    i = res - 1
    return arr[i] if 0 <= i < len(arr) else None


def win(arr, res, half):
    if arr is None:
        return None
    seg = arr[max(0, res - 1 - half) : min(len(arr), res + half)]
    return sum(seg) / len(seg) if seg else None


def coil_run(ss3, res):
    if not ss3:
        return None
    i = res - 1
    if not (0 <= i < len(ss3)):
        return None
    if ss3[i] != "C":
        return 0
    lo = hi = i
    while lo > 0 and ss3[lo - 1] == "C":
        lo -= 1
    while hi + 1 < len(ss3) and ss3[hi + 1] == "C":
        hi += 1
    return hi - lo + 1


def segments(prt: str, ch: str):
    out, start = [], None
    for i, c in enumerate(prt):
        if c == ch:
            if start is None:
                start = i
        elif start is not None:
            out.append((start + 1, i))
            start = None
    if start is not None:
        out.append((start + 1, len(prt)))
    return out


def translate_orf(cds: str | None) -> str | None:
    """L1 ORF -> protein. Translation stops at the first stop codon, which is
    also how a trailing stop is dropped; a CDS whose length is not a multiple
    of three is truncated to the last whole codon rather than silently
    frame-shifting."""
    if not cds:
        return None
    seq = cds.strip().upper()
    seq = seq[: len(seq) - len(seq) % 3]
    if not seq:
        return None
    return str(Seq(seq).translate(to_stop=True))


def compare_sequences(l1: str | None, stored: str) -> tuple[str, str]:
    """Exact comparison, because equal length is not the same protein.

    Returns a verdict and a human-readable detail. Anything other than
    ``exact`` means the stored per-residue features describe a different
    sequence from the one L1 builds, so its residue numbers do not transfer.
    """
    if not stored:
        return "no stored sequence", ""
    if l1 is None:
        return "no L1 ORF", ""
    if l1 == stored:
        return "exact", ""
    if len(l1) != len(stored):
        return "length differs", f"L1 {len(l1)} aa vs stored {len(stored)} aa"
    first = next(i for i in range(len(l1)) if l1[i] != stored[i])
    n_diff = sum(1 for a, b in zip(l1, stored, strict=True) if a != b)
    return (
        "same length, different sequence",
        f"{n_diff} substitution(s), first at residue {first + 1}: "
        f"L1 {l1[first]} vs stored {stored[first]}",
    )


def git_sha() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=HERE, text=True).strip()
    except Exception:
        return "unknown"


def parse_manifest(path: Path):
    """candidate rows -> one record per distinct proteoform."""
    rows = [ln.rstrip("\n").split("\t") for ln in path.open()]
    h = {c: i for i, c in enumerate(rows[0])}
    by_key: dict[tuple, dict] = {}
    for r in rows[1:]:
        key, kind_set, species = r[h["key"]], r[h["set"]], r[h["species"]]
        parts = key.split(":")
        if parts[0] == "lib":
            ident, gene, acc = ("gene", parts[1]), parts[1], None
        elif parts[0] == "iso":
            ident, gene, acc = ("acc", parts[2]), parts[1], parts[2]
        elif parts[0] == "ort":
            ident, gene, acc = ("acc", parts[3]), parts[2], parts[3]
        else:
            continue
        rec = by_key.setdefault(
            ident,
            {
                "key": key,
                "gene": gene,
                "acc": acc,
                "set": kind_set,
                "species": species,
                "transcript": r[h["transcript"]],
                "candidates": [],
            },
        )
        rec["candidates"].append(r[h["candidate_id"]])
    return list(by_key.values())


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--manifest", required=True, type=Path)
    ap.add_argument(
        "--cds", required=True, type=Path, help="candidate_cds.json — the L1 ORFs themselves"
    )
    ap.add_argument("--features", default="site_features.tsv")
    ap.add_argument("--cache", type=Path, default=None, help="feature cache JSONL")
    ap.add_argument("--env-from", default=None)
    ap.add_argument("--terminal-margin", type=int, default=25)
    ap.add_argument("--n-sites", type=int, default=3)
    ap.add_argument("--min-sep", type=int, default=25)
    ap.add_argument("--batch", type=int, default=25)
    ap.add_argument("--limit", type=int, default=0, help="score only the first N proteoforms")
    ap.add_argument("--seed", type=int, default=20261014)
    ap.add_argument("--write-d1", action="store_true", help="publish to tag_site_prediction")
    ap.add_argument("--run-id", default="l1_library_v1")
    ap.add_argument("--out", default="l1_library_sites.tsv")
    args = ap.parse_args()
    np.random.seed(args.seed)

    try:
        from deliverome_analysis.env import load_env

        load_env()
    except Exception:
        pass
    if args.env_from:
        from dotenv import load_dotenv

        load_dotenv(Path(args.env_from).expanduser(), override=False)

    # The L1 ORFs themselves. Joined on (gene, transcript) because the CDS
    # file keys isoforms by transcript while the manifest keys them by
    # accession.
    cds_by_tx: dict[tuple[str, str], str] = {}
    for k, v in json.loads(args.cds.read_text()).items():
        parts = k.split(":")
        gene = parts[1] if parts[0] != "ort" else parts[2]
        cds_by_tx[(gene, v["transcript"])] = v["cds"]
    print(f"L1 ORFs loaded: {len(cds_by_tx)}")

    proteoforms = parse_manifest(args.manifest)
    if args.limit:
        proteoforms = proteoforms[: args.limit]
    print(
        f"L1: {sum(len(p['candidates']) for p in proteoforms)} candidates, "
        f"{len(proteoforms)} distinct proteoforms"
    )

    # lib:GENE rows name no accession; resolve each to its canonical one.
    genes = sorted({p["gene"] for p in proteoforms if p["acc"] is None})
    gene_acc: dict[str, str] = {}
    for i in range(0, len(genes), 60):
        b = genes[i : i + 60]
        ph = ",".join("?" * len(b))
        for r in d1(
            "SELECT DISTINCT gene_symbol, uniprot_acc FROM topology_public "
            f"WHERE tool_version LIKE 'deeptmhmm-1%' AND gene_symbol IN ({ph});",
            b,
        ):
            gene_acc.setdefault(r["gene_symbol"], r["uniprot_acc"])
    for p in proteoforms:
        if p["acc"] is None:
            p["acc"] = gene_acc.get(p["gene"])
        p["is_canonical"] = "-" not in (p["acc"] or "") and p["set"] in (
            "library_reference",
            "expression_rescue",
            "giant_truncation",
        )
    print(f"resolved {len(gene_acc)}/{len(genes)} library_reference genes to an accession")

    accs = sorted({p["acc"] for p in proteoforms if p["acc"]})

    cache: dict[str, dict] = {}
    cache_path = args.cache or (HERE / "l1_feature_cache.jsonl")
    if cache_path.exists():
        for ln in cache_path.open():
            rec = json.loads(ln)
            cache[rec["acc"]] = rec
        print(f"cache: {len(cache)} proteoforms already fetched")

    todo = [a for a in accs if a not in cache]
    print(f"fetching topology + disorder for {len(todo)} proteoforms")
    with cache_path.open("a") as fh:
        for i in range(0, len(todo), args.batch):
            b = todo[i : i + args.batch]
            ph = ",".join("?" * len(b))
            topo = {
                r["uniprot_acc_full"]: r
                for r in d1(
                    "SELECT uniprot_acc_full, gene_symbol, protein_length, sequence, "
                    "per_residue_topology, tm_helix_count, signal_peptide_length, "
                    "ecd_length_residues, icd_length_residues FROM topology_public "
                    f"WHERE tool_version LIKE 'deeptmhmm-1%' AND uniprot_acc_full IN ({ph});",
                    b,
                )
            }
            dis: dict[str, dict] = defaultdict(dict)
            for r in d1(
                "SELECT uniprot_acc_full, predictor, disorder_hex, rsa_hex, plddt_hex, ss3 "
                f"FROM disorder_public WHERE uniprot_acc_full IN ({ph});",
                b,
            ):
                dis[r["uniprot_acc_full"]][r["predictor"]] = r
            for a in b:
                rec = {"acc": a, "topo": topo.get(a), "dis": dis.get(a, {})}
                cache[a] = rec
                fh.write(json.dumps(rec) + "\n")
            if (i // args.batch) % 10 == 0:
                print(f"  {min(i + args.batch, len(todo))}/{len(todo)}", flush=True)

    # Published predictors, canonical coordinates only. Small table; take it whole.
    pub: dict[tuple[str, str], dict] = {}
    for pred in ("prodomino", "epictope", "indeli_e"):
        for r in d1(
            "SELECT uniprot_acc, protein_length, scores_hex, score_min, score_max "
            "FROM tag_site_prediction WHERE predictor = ?;",
            [pred],
        ):
            pub[(pred, r["uniprot_acc"])] = r
    print(f"published predictor arrays available for {len(pub)} (predictor, protein) pairs")

    matched = pd.read_csv(HERE / args.features, sep="\t", low_memory=False)
    for c in ("ss3", "topology_class", "region_class"):
        matched[c] = matched[c].fillna("").astype(str)
    y = matched["label"].to_numpy()
    model = models(int(y.sum()), FUSED_NUMERIC)["logistic"]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model.fit(matched, y)

    out_rows, arrays, skipped = [], [], []
    for n, p in enumerate(proteoforms):
        acc = p["acc"]
        rec = cache.get(acc) if acc else None
        if not rec or not rec.get("topo"):
            skipped.append(
                {
                    **{k: p[k] for k in ("key", "gene", "set", "species")},
                    "acc": acc,
                    "reason": "no stored topology",
                }
            )
            continue
        t = rec["topo"]
        prt = t.get("per_residue_topology") or ""
        plen = int(t.get("protein_length") or len(prt) or 0)
        l1_protein = translate_orf(cds_by_tx.get((p["gene"], p["transcript"])))
        stored_seq = (t.get("sequence") or "").strip().upper()
        verdict, detail = compare_sequences(l1_protein, stored_seq)
        if verdict != "exact":
            skipped.append(
                {
                    **{k: p[k] for k in ("key", "gene", "set", "species")},
                    "acc": acc,
                    "reason": verdict,
                    "detail": detail,
                }
            )
            continue
        s = rec["dis"]
        meta = s.get("metapredict") or {}
        nsp = s.get("netsurfp-3.0") or {}
        afd = s.get("alphafold-disorder") or {}
        dis = unit(meta.get("disorder_hex"))
        rsa = unit(afd.get("rsa_hex")) or unit(nsp.get("rsa_hex"))
        pl = plddt(afd.get("plddt_hex"))
        ss3 = nsp.get("ss3")
        if dis is None and rsa is None and pl is None:
            skipped.append(
                {
                    **{k: p[k] for k in ("key", "gene", "set", "species")},
                    "acc": acc,
                    "reason": "no disorder/RSA/pLDDT",
                }
            )
            continue

        pub_arrays = {}
        if p["is_canonical"]:
            for pred in ("prodomino", "epictope", "indeli_e"):
                r = pub.get((pred, acc))
                if r and int(r["protein_length"] or 0) == plen:
                    pub_arrays[pred] = [b / 255.0 for b in bytes.fromhex(r["scores_hex"])]

        tm = segments(prt, "M")
        osegs = segments(prt, "O")
        recs = []
        for pos in range(1, plen + 1):
            if not (0 <= pos - 1 < len(prt)) or prt[pos - 1] != "O":
                continue
            if min(pos, plen - pos + 1) <= args.terminal_margin:
                continue
            d_at = at(dis, pos)
            row = {c: np.nan for c in CANONICAL_ONLY}
            for pred, arr in pub_arrays.items():
                row[pred] = at(arr, pos)
                row[f"{pred}_w5"] = win(arr, pos, 5)
            row.update(
                {
                    "disorder": d_at,
                    "disorder_w5": win(dis, pos, 5),
                    "disorder_w15": win(dis, pos, 15),
                    "rsa": at(rsa, pos),
                    "rsa_w5": win(rsa, pos, 5),
                    "plddt": at(pl, pos),
                    "plddt_w5": win(pl, pos, 5),
                    "ss3": ss3[pos - 1] if ss3 and 0 <= pos - 1 < len(ss3) else "",
                    "loop_len": coil_run(ss3, pos),
                    "tm_helix_count": t.get("tm_helix_count"),
                    "signal_peptide_length": t.get("signal_peptide_length"),
                    "ecd_length_residues": t.get("ecd_length_residues"),
                    "icd_length_residues": t.get("icd_length_residues"),
                    "dist_to_tm": (
                        min(0 if a <= pos <= b else min(abs(pos - a), abs(pos - b)) for a, b in tm)
                        if tm
                        else None
                    ),
                    "dist_to_terminus": min(pos, plen - pos + 1),
                    "rel_position": pos / plen,
                    "protein_length": plen,
                    "topology_class": "multipass_both_cyto",
                    "region_class": (
                        "" if d_at is None else ("disordered" if d_at >= 0.5 else "ordered_surface")
                    ),
                    "residue": pos,
                    "region_index": next(
                        (i for i, (a, b) in enumerate(osegs) if a <= pos <= b), None
                    ),
                }
            )
            recs.append(row)
        if not recs:
            skipped.append(
                {
                    **{k: p[k] for k in ("key", "gene", "set", "species")},
                    "acc": acc,
                    "reason": "no scorable extracellular residue",
                }
            )
            continue

        df = pd.DataFrame(recs)
        for c in FUSED_NUMERIC:
            if c not in df.columns:
                df[c] = np.nan
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            df["score"] = model.predict_proba(df)[:, 1]

        full = np.zeros(plen, dtype=float)
        full[df["residue"].to_numpy() - 1] = df["score"].to_numpy()
        arrays.append(
            {
                "acc": acc,
                "gene": t.get("gene_symbol") or p["gene"],
                "length": plen,
                "scores": full,
                "n_scored": len(df),
                "predictors": ";".join(sorted(pub_arrays)) or "",
            }
        )

        df = df.sort_values("score", ascending=False)
        picked = []
        for _, row in df.iterrows():
            if len(picked) >= args.n_sites:
                break
            if any(abs(row["residue"] - q["residue"]) < args.min_sep for q in picked):
                continue
            picked.append(row)
        for rank, row in enumerate(picked, 1):
            out_rows.append(
                {
                    "candidate_ids": ";".join(p["candidates"]),
                    "key": p["key"],
                    "gene": p["gene"],
                    "uniprot_acc_full": acc,
                    "set": p["set"],
                    "species": p["species"],
                    "protein_length": plen,
                    "rank": rank,
                    "insert_after_residue": int(row["residue"]),
                    "score": round(float(row["score"]), 4),
                    "within_protein_percentile": round(
                        float((df["score"] < row["score"]).mean()), 4
                    ),
                    "loop_len": row.get("loop_len"),
                    "rsa": row.get("rsa"),
                    "disorder": row.get("disorder"),
                    "plddt": row.get("plddt"),
                    "n_extracellular_scored": len(df),
                    "published_predictors": ";".join(sorted(pub_arrays)) or "none",
                }
            )
        if (n + 1) % 500 == 0:
            print(f"  scored {n + 1}/{len(proteoforms)}", flush=True)

    out = pd.DataFrame(out_rows)
    out.to_csv(HERE / args.out, sep="\t", index=False)
    pd.DataFrame(skipped).to_csv(HERE / "l1_library_skipped.tsv", sep="\t", index=False)
    print(f"\nscored {len(arrays)} proteoforms, {len(out)} sites -> {args.out}")
    print(f"skipped {len(skipped)} -> l1_library_skipped.tsv")
    if skipped:
        print(
            pd.DataFrame(skipped)["reason"]
            .str.replace(r"\d+", "N", regex=True)
            .value_counts()
            .to_string()
        )

    prov = {
        "run_id": args.run_id,
        "seed": args.seed,
        "git_commit": git_sha(),
        "generated_at": datetime.now(UTC).isoformat(),
        "manifest": str(args.manifest),
        "proteoforms_scored": len(arrays),
        "proteoforms_skipped": len(skipped),
        "predictor": PREDICTOR,
        "predictor_version": PREDICTOR_VERSION,
        "with_published_predictors": sum(1 for a in arrays if a["predictors"]),
    }
    (HERE / "l1_library_provenance.json").write_text(json.dumps(prov, indent=2) + "\n")

    if args.write_d1:
        print(
            f"\npublishing {len(arrays)} rows to tag_site_prediction as {PREDICTOR} "
            f"v{PREDICTOR_VERSION}, run_id {args.run_id}"
        )
        d1(
            "DELETE FROM tag_site_prediction WHERE predictor = ? AND run_id = ?;",
            [PREDICTOR, args.run_id],
        )
        now = datetime.now(UTC).isoformat()
        cfg = (
            "fused logistic over topology/disorder/RSA/pLDDT/loop features; "
            "published predictors used only where the proteoform is canonical"
        )
        sent = 0
        for i in range(0, len(arrays), 20):
            chunk = arrays[i : i + 20]
            sql = (
                "INSERT INTO tag_site_prediction (uniprot_acc, gene_symbol, predictor, "
                "predictor_version, model_config, protein_length, n_scores, scores_hex, "
                "score_min, score_max, licence, citation, run_id, created_at) VALUES "
                + ",".join(["(?,?,?,?,?,?,?,?,?,?,?,?,?,?)"] * len(chunk))
                + ";"
            )
            params = []
            for a in chunk:
                sc = a["scores"]
                params += [
                    a["acc"],
                    a["gene"],
                    PREDICTOR,
                    PREDICTOR_VERSION,
                    cfg,
                    a["length"],
                    int(a["n_scored"]),
                    bytes(round(min(max(v, 0.0), 1.0) * 255) for v in sc).hex(),
                    float(sc.min()),
                    float(sc.max()),
                    "CC-BY-4.0",
                    "Deliverome S5E13 tag-site model",
                    args.run_id,
                    now,
                ]
            d1(sql, params)
            sent += len(chunk)
            if (i // 20) % 25 == 0:
                print(f"  {sent}/{len(arrays)}", flush=True)
        print(f"published {sent} rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
