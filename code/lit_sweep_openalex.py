#!/usr/bin/env python3
"""WP 1.2 supplement: OpenAlex forward-citation sweep (rate-limit friendly).

Covers both index records of the two closest papers (arXiv version + published
version) and writes screened shortlists to data/lit/.
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

MAILTO = "scf-sweep@example.org"

TARGETS = {
    # label: list of OpenAlex work ids whose citing works we pull
    "cevid": ["W2900832492", "W3117548431"],
    "rendsburg": ["W2211.01903_PLACEHOLDER"],
}

SCREEN_RE = re.compile(
    r"(phase transition|detectab|minimax|lower bound|le cam|adaptive tun|"
    r"empirical bayes|frontier|tracy.?widom|bbp|baik|benaych|johnstone|"
    r"bias of|debias|outlier eigenvalue|spiked|spike model|factor structure|"
    r"confounding strength)", re.I)


def oa_get(url):
    sep = "&" if "?" in url else "?"
    url = f"{url}{sep}mailto={MAILTO}"
    for attempt in range(5):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": f"scf-phase1/0.1 ({MAILTO})"})
            with urllib.request.urlopen(req, timeout=30) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            print(f"[oa] {e}; retry", flush=True)
            time.sleep(5 * (attempt + 1))
    return None


def find_id_by_title(title):
    q = urllib.parse.quote(title)
    d = oa_get(f"https://api.openalex.org/works?filter=title.search:{q}&per-page=5")
    if not d:
        return None
    best = d.get("results", [])
    return best[0]["id"].rsplit("/", 1)[-1] if best else None


def pull_citations(label, work_ids):
    all_rows = []
    for wid in work_ids:
        cursor = "*"
        while cursor:
            d = oa_get(f"https://api.openalex.org/works?filter=cites:{wid}"
                       f"&per-page=200&cursor={cursor}"
                       "&select=id,title,publication_year,cited_by_count,doi,"
                       "primary_location,abstract_inverted_index")
            if not d:
                break
            for w in d.get("results", []):
                abstract = ""
                aii = w.get("abstract_inverted_index") or {}
                if aii:
                    pos = {}
                    for word, idxs in aii.items():
                        for i in idxs:
                            pos[i] = word
                    abstract = " ".join(pos[i] for i in sorted(pos))
                blob = f"{w.get('title') or ''} {abstract}"
                src = ((w.get("primary_location") or {}).get("source") or {})
                all_rows.append({
                    "from_work": wid,
                    "year": w.get("publication_year"),
                    "title": w.get("title") or "",
                    "venue": src.get("display_name") or "",
                    "doi": w.get("doi") or "",
                    "cites": w.get("cited_by_count"),
                    "screen_hit": bool(SCREEN_RE.search(blob)),
                    "hit_terms": ",".join(sorted({m.group(0).lower()
                                                  for m in SCREEN_RE.finditer(blob)})),
                })
            cursor = (d.get("meta") or {}).get("next_cursor")
            time.sleep(0.5)
    out = os.path.join(BASE, f"openalex_fwd_{label}.json")
    seen = set()
    uniq = []
    for r_ in all_rows:
        key = r_["title"].lower(), r_["year"]
        if key in seen:
            continue
        seen.add(key)
        uniq.append(r_)
    with open(out, "w") as f:
        json.dump(uniq, f, indent=1)
    print(f"[{label}] {len(uniq)} unique citing works "
          f"({sum(r_['screen_hit'] for r_ in uniq)} keyword hits)", flush=True)


if __name__ == "__main__":
    cevid_ids = TARGETS["cevid"]
    rid = find_id_by_title("A consistent estimator for confounding strength")
    print("[rendsburg] openalex id:", rid, flush=True)
    pull_citations("cevid", [i for i in cevid_ids if i])
    if rid:
        pull_citations("rendsburg", [rid])
