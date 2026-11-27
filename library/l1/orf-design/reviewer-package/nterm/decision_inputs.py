"""Every S3 input the N-after-SP decisions and arms are built from, pinned by SHA-256.

``nterm_decisions.tsv`` and ``constructs_l1_nterm.*`` are deterministic functions of these
files and the committed code: ``nterm_decide.py decide`` then ``build_l1_nterm_arms.py --check``
reproduce them byte for byte. Two inputs are republished by other sessions, so this script
freezes a copy of each under a content-addressed key of ours before recording it:

- l1-tag-refold's ``verdicts_nterm.tsv`` (``refold-verdicts/verdicts_nterm.<sha12>.tsv``);
- nothing else today; every SignalP run, the D1 snapshot, the ortholog rows and the UniProt
  pro-peptide snapshot are written once under dated or run-named keys.

Tool versions are in the files themselves: SignalP 6.0 package and model fingerprints in each
run's ``provenance.json`` (``signalp6_identity.json`` pins them), disorder predictor and
DeepTMHMM v1 versions in the D1 rows, UniProt release in ``uniprot_features_snapshot.json``,
ColabFold settings in l1-tag-refold's ``structures.tsv``.

Written by the Claude Code session nterm-reviewer-package (PR #254), 2026-11-13, for the
overnight design review relayed from Becca.

    uv run --locked python decision_inputs.py            # freeze and write the manifest
    uv run --locked python decision_inputs.py --verify   # rehash every pinned key
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))

import nterm_decide as nd  # noqa: E402
import nterm_overnight as no  # noqa: E402
import nterm_placement as npl  # noqa: E402
import propeptide_check as ppc  # noqa: E402

from deliverome_analysis.cloning import library_design as ld  # noqa: E402

OUT = HERE / "decision_inputs.json"
FROZEN = npl.PREFIX + "reviewer-package/nterm/refold-verdicts/"


def _sha(path: str) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _get(key: str) -> str:
    from deliverome_analysis.core.s3 import external, processed

    if key.startswith("external://"):
        return external.download(key.removeprefix("external://"))
    return processed.download(key)


def keys() -> dict[str, str]:
    k = {
        "designed CDS": ld.CDS_KEY,
        "domesticated ORFs": ld.CONSTRUCTS_KEY,
        "D1 predictor rows (SignalP, disorder, input sequences)": npl.ROWS_KEY,
        "interval census": npl.CENSUS_KEY,
        "re-prediction sweep submission": npl.SWEEP_SUBMISSION_KEY,
        "re-prediction sweep calls": npl.SWEEP_REPREDICTION_KEY,
        "ortholog predictor rows": no.ORTHO_KEY,
        "UniProt pro-peptide snapshot": json.loads(
            (HERE / "uniprot_features_snapshot.json").read_text()
        )["key"],
    }
    sp = HERE / "uniprot_sp_disagree_snapshot.json"
    if sp.exists():
        k["UniProt second-source snapshot (signal peptide, TM, topology)"] = json.loads(
            sp.read_text()
        )["key"]
    for run in no.RUN_IDS:
        k[f"SignalP run {run} calls"] = no.calls_key(run)
        k[f"SignalP run {run} provenance"] = no.calls_key(run).replace(
            "calls.json", "provenance.json"
        )
    return k


def main(verify: bool) -> int:
    from deliverome_analysis.core.s3 import processed

    if verify:
        pinned = json.loads(OUT.read_text())["inputs"]
        bad = [
            n
            for n, e in pinned.items()
            if _sha(str(HERE / e["key"][4:]) if e["key"].startswith("git:") else _get(e["key"]))
            != e["sha256"]
        ]
        print(f"{len(pinned) - len(bad)}/{len(pinned)} inputs match their pinned SHA-256")
        for n in bad:
            print("  changed:", n)
        return 1 if bad else 0
    inputs = {}
    for name, key in keys().items():
        if not key.startswith("external://") and not processed.exists(key):
            continue
        inputs[name] = {"key": key, "sha256": _sha(_get(key))}
    used = json.loads((HERE / "nterm_decisions_summary.json").read_text())
    source = used.get("refold_verdicts_source", nd.VERDICTS_KEY)
    if (HERE / source).exists():  # a table committed on main (tag_refold_verdicts.tsv)
        assert _sha(str(HERE / source)) == used["refold_verdicts_sha256"]
        inputs["refold verdicts (committed by l1-tag-refold, PR #301)"] = {
            "key": "git:" + source,
            "sha256": used["refold_verdicts_sha256"],
        }
    else:
        v = processed.download(source)
        vsha = _sha(v)
        frozen = FROZEN + f"verdicts_nterm.{vsha[:12]}.tsv"
        if not processed.exists(frozen):
            processed.upload(v, frozen)
        inputs["refold verdicts (frozen copy of l1-tag-refold's verdicts_nterm.tsv)"] = {
            "key": frozen,
            "sha256": vsha,
            "source": source,
        }
    manifest = {
        "inputs": inputs,
        "uniprot_release": json.loads((HERE / "uniprot_features_snapshot.json").read_text())[
            "release"
        ],
        "signalp_identity": "signalp6_identity.json",
        "lddt_borderline_margin": nd.LDDT_MARGIN,
        "activation_only_propeptides": sorted(ppc.ACTIVATION_ONLY),
    }
    OUT.write_text(json.dumps(manifest, indent=1) + "\n")
    print(f"{len(inputs)} inputs pinned -> {OUT.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main("--verify" in sys.argv))
