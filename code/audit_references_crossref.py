"""Audit all manuscript references against Crossref metadata without rewriting citations automatically."""
from __future__ import annotations

import json
import re
import time
import os
from pathlib import Path
from urllib.parse import quote
from urllib.request import Request, urlopen

import pandas as pd
from docx import Document

ROOT = Path(__file__).resolve().parents[1]
SOURCE = Path(os.environ["HVI_MANUSCRIPT_SOURCE"]).expanduser()
OUT = ROOT / "results"


def clean(s: str) -> str:
    return re.sub(r"[^a-z0-9 ]+", " ", s.lower()).strip()


def title_tokens(s: str) -> set[str]:
    return {x for x in clean(s).split() if len(x) > 3}


def extract_refs():
    d = Document(SOURCE)
    refs=[]
    for p in d.paragraphs[158:]:
        text=p.text.strip().replace("\n"," ")
        if not text or text.startswith("Supporting information"):
            break
        refs.append(text)
    return refs


def get_json(url):
    req=Request(url, headers={"User-Agent":"Kenya-HVI-resubmission-audit/1.0 (reference verification)"})
    with urlopen(req, timeout=25) as r:
        return json.loads(r.read().decode("utf-8"))


def year_from(text):
    m=re.search(r"\b(19|20)\d{2}\b", text)
    return int(m.group()) if m else None


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    rows=[]
    for num, ref in enumerate(extract_refs(), start=1):
        doi_match=re.search(r"doi:\s*(10\.\S+)", ref, flags=re.I)
        doi=doi_match.group(1).rstrip(".") if doi_match else ""
        try:
            if doi:
                item=get_json("https://api.crossref.org/works/"+quote(doi,safe=""))["message"]
                retrieval="DOI lookup"
            else:
                data=get_json("https://api.crossref.org/works?rows=1&query.bibliographic="+quote(ref))
                items=data["message"].get("items",[])
                item=items[0] if items else {}
                retrieval="bibliographic search"
            matched_title=(item.get("title") or [""])[0]
            matched_doi=item.get("DOI","")
            matched_year=((item.get("published-print") or item.get("published-online") or item.get("issued") or {}).get("date-parts") or [[None]])[0][0]
            submitted_year=year_from(ref)
            overlap=len(title_tokens(ref)&title_tokens(matched_title))/max(1,len(title_tokens(matched_title)))
            if doi and matched_doi.lower()==doi.lower(): status="verified DOI"
            elif overlap>=.75 and (submitted_year is None or matched_year is None or abs(submitted_year-matched_year)<=1): status="probable match; manual check"
            elif item: status="possible mismatch; manual check"
            else: status="not found in Crossref; manual publisher check"
            rows.append({"reference_number":num,"submitted_reference":ref,"submitted_year":submitted_year,"submitted_doi":doi,"retrieval":retrieval,"crossref_status":status,"crossref_title":matched_title,"crossref_year":matched_year,"crossref_doi":matched_doi,"title_token_overlap":round(overlap,3),"crossref_score":item.get("score","")})
        except Exception as e:
            rows.append({"reference_number":num,"submitted_reference":ref,"submitted_year":year_from(ref),"submitted_doi":doi,"retrieval":"error","crossref_status":"manual check required (lookup error)","crossref_title":"","crossref_year":"","crossref_doi":"","title_token_overlap":"","crossref_score":"","error":str(e)[:180]})
        time.sleep(.12)
    audit=pd.DataFrame(rows)
    audit.to_csv(OUT/"reference_verification_audit_crossref.csv",index=False)
    print(audit.crossref_status.value_counts().to_string())
    print("references audited",len(audit))

if __name__=="__main__":
    main()
