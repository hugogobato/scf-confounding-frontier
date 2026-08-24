#!/usr/bin/env python3
"""WP 1.2 literature sweep: forward citations + locator pinning + arXiv queries.

Pulls raw API responses into data/lit/ so the novelty memo is source-traceable.
Paced for unauthenticated Semantic Scholar limits; retries on 429.
"""
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

BASE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "lit")
os.makedirs(BASE, exist_ok=True)

S2 = "https://api.semanticscholar.org/graph/v1"
FIELDS = "title,year,venue,externalIds,citationCount,abstract"


def s2_get(url, tag, max_tries=12):
    """GET with 429 backoff; cache raw JSON to disk."""
    path = os.path.join(BASE, tag + ".json")
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    for attempt in range(max_tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "scf-phase1/0.1"})
            with urllib.request.urlopen(req, timeout=30) as r:
                data = json.loads(r.read().decode())
            with open(path, "w") as f:
                json.dump(data, f)
            time.sleep(4)
            return data
        except Exception as e:
            wait = min(30, 5 * (attempt + 1))
            print(f"[{tag}] {e}; retry in {wait}s", flush=True)
            time.sleep(wait)
    raise RuntimeError(f"failed: {url}")


CITED = {
    "cevid_1811.05352": "arXiv:1811.05352",
    "rendsburg_2211.01903": "arXiv:2211.01903",
}


def pull_forward_citations():
    for name, pid in CITED.items():
        offset = 0
        limit = 100
        while True:
            url = (f"{S2}/paper/{pid}/citations?fields={FIELDS}"
                   f"&limit={limit}&offset={offset}")
            data = s2_get(url, f"fwd_{name}_off{offset}")
            citing = data.get("data") or []
            n_next = data.get("next")
            print(f"[{name}] offset={offset} got={len(citing)} next={n_next}", flush=True)
            if not citing or n_next is None or n_next <= offset:
                break
            offset = n_next


ARXIV_QUERIES = [
    'all:"deconfounding" AND all:"phase transition"',
    'all:"hidden confounding" AND all:"random matrix"',
    'all:"deconfounding" AND all:"bias"',
    'all:"spectral deconfounding"',
    'all:"confounding" AND all:"Benaych-Georges"',
    'all:"confounding" AND all:"Tracy-Widom"',
    'all:"confounder detection" AND all:"spike"',
    'all:"detectability frontier"',
    'all:"confounding strength" AND all:"random matrix"',
]


def pull_arxiv_queries():
    for i, q in enumerate(ARXIV_QUERIES):
        path = os.path.join(BASE, f"arxiv_q{i}.txt")
        if os.path.exists(path):
            continue
        url = ("http://export.arxiv.org/api/query?search_query="
               + urllib.parse.quote(q)
               + "&start=0&max_results=50&sortBy=submittedDate&sortOrder=descending")
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "scf-phase1/0.1"})
            with urllib.request.urlopen(req, timeout=30) as r:
                text = r.read().decode()
            with open(path, "w") as f:
                f.write(f"# QUERY: {q}\n" + text)
            print(f"[arxiv q{i}] saved ({len(text)} bytes): {q}", flush=True)
        except Exception as e:
            print(f"[arxiv q{i}] FAILED: {e}", flush=True)
        time.sleep(3)


# Classical sources whose locators must be pinned (WP 1.2 action 2).
LOCATOR_TARGETS = [
    ("bbp_2005", "Phase transition in the largest eigenvalue of covariance matrices"),
    ("johnstone_2001", "On the distribution of the largest eigenvalue in principal components analysis"),
    ("bgn_2011", "The eigenvalues and eigenvectors of finite, low rank perturbations of large random matrices"),
    ("bgn_2012", "The singular values and vectors of low rank perturbations of large rectangular random matrices"),
    ("omh_2013", "Asymptotic power of sphericity tests for high-dimensional data"),
    ("onatski_2010", "Determining the number of factors from empirical distribution of eigenvalues"),
    ("fan_liao_2014", "Endogeneity in high dimensions"),
    ("wang_blei_2019", "The blessings of multiple causes"),
    ("bai_ng_2002", "Determining the number of factors in approximate factor models"),
    ("dobriban_wager_2018", "High-dimensional asymptotics of prediction: Ridge regression and classification"),
    ("hastie_2022", "Surprises in high-dimensional ridgeless least squares interpolation"),
    ("knowles_yin_2017", "Anisotropic local laws for random matrices"),
    ("janzing_scholkopf_2018", "Detecting confounding in multivariate linear models via spectral analysis"),
    ("cevid_2020_aos", "Spectral Deconfounding via Perturbed Sparse Linear Models"),
    ("rendsburg_2022", "A Consistent Estimator for Confounding Strength"),
    ("schur_peters_2024", "Deconfounding time series with robust regression"),
    ("scheidegger_2024", "Spectral deconfounding for high-dimensional sparse additive models"),
    ("ulmer_2025", "Spectrally deconfounded random forests"),
    ("chernozhukov_lava", "High-dimensional sparse econometric models LAVA"),
]


def pin_locators():
    """Pin DOIs via Crossref (reliable for classical journal sources)."""
    results = {}
    path = os.path.join(BASE, "locator_pins.json")
    if os.path.exists(path):
        return
    for key, title in LOCATOR_TARGETS:
        url = ("https://api.crossref.org/works?query.bibliographic="
               + urllib.parse.quote(title) + "&rows=3&mailto=scf-sweep@example.org")
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "scf-phase1/0.1 (mailto:scf-sweep@example.org)"})
            with urllib.request.urlopen(req, timeout=30) as r:
                data = json.loads(r.read().decode())
            items = data.get("message", {}).get("items", [])
            best = []
            for it in items[:3]:
                best.append({
                    "title": (it.get("title") or [""])[0],
                    "DOI": it.get("DOI"),
                    "container": (it.get("container-title") or [""])[0],
                    "year": (it.get("issued", {}).get("date-parts") or [[None]])[0][0],
                    "volume": it.get("volume"),
                    "issue": it.get("issue"),
                    "page": it.get("page"),
                })
            results[key] = {"query": title, "candidates": best}
            if best:
                b = best[0]
                print(f"[pin:{key}] -> {b['title']!r} doi:{b['DOI']} {b['container']} "
                      f"{b['volume']}({b['issue']}) {b['year']}", flush=True)
            else:
                print(f"[pin:{key}] NO CROSSREF RESULT", flush=True)
        except Exception as e:
            print(f"[pin:{key}] FAILED: {e}", flush=True)
        time.sleep(1)
    with open(path, "w") as f:
        json.dump(results, f, indent=1)


SCREEN_RE = re.compile(
    r"(phase transition|detectab|minimax|lower bound|le cam|adaptive tun|"
    r"empirical bayes|frontier|tracy.?widom|bbp|baik|benaych|johnstone|"
    r"bias of|debias|outlier eigenvalue|spiked|spike model|factor structure)", re.I)


def screen():
    rows = []
    for name in CITED:
        offset = 0
        while True:
            path = os.path.join(BASE, f"fwd_{name}_off{offset}.json")
            if not os.path.exists(path):
                break
            with open(path) as f:
                data = json.load(f)
            for rec in data.get("data") or []:
                p = rec.get("citingPaper") or {}
                blob = " ".join(str(p.get(k) or "") for k in ("title", "abstract"))
                ext = p.get("externalIds") or {}
                rows.append({
                    "source": name,
                    "year": p.get("year"),
                    "title": (p.get("title") or "").strip(),
                    "venue": p.get("venue") or "",
                    "doi": ext.get("DOI", ""),
                    "arxiv": ext.get("ArXiv", ""),
                    "cites": p.get("citationCount"),
                    "screen_hit": bool(SCREEN_RE.search(blob)),
                    "hit_terms": ",".join(sorted(set(m.group(0).lower() for m in SCREEN_RE.finditer(blob)))),
                })
            nxt_path = os.path.join(BASE, f"fwd_{name}_off{offset + 100}.json")
            if not os.path.exists(nxt_path):
                break
            offset += 100
    out = os.path.join(BASE, "forward_citations_screened.json")
    with open(out, "w") as f:
        json.dump(rows, f, indent=1)
    hits = [r for r in rows if r["screen_hit"]]
    hits.sort(key=lambda r: -(r["cites"] or 0))
    with open(os.path.join(BASE, "shortlist.txt"), "w") as f:
        for r in hits:
            f.write(f"{r['year']} | cites={r['cites']} | {r['title']} | "
                    f"arXiv:{r['arxiv'] or '-'} | doi:{r['doi'] or '-'} | [{r['hit_terms']}]\n")
    print(f"screened {len(rows)} citing papers; {len(hits)} keyword hits -> shortlist.txt", flush=True)


if __name__ == "__main__":
    step = sys.argv[1] if len(sys.argv) > 1 else "all"
    if step in ("all", "pins"):
        pin_locators()
    if step in ("all", "fwd"):
        pull_forward_citations()
    if step in ("all", "arxiv"):
        pull_arxiv_queries()
    if step == "screen":
        screen()
    print("done:", step, flush=True)
