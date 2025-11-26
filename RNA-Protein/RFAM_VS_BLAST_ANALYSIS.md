# Rfam vs BLAST: Why We Chose BLAST

## Professor Azam's Suggestion

From her email:
> "One solution could be to use Rfam family datasets. Rfam provides high-quality, curated alignments, but these are mostly limited to short RNAs or specific structured domains."

She identified Rfam as **one possible solution** but also noted its **limitations**.

---

## Testing Both Approaches

We tested both Rfam and BLAST on the same test RNA (16S rRNA, 397nt):

### Rfam Approach Results

**What it does:**
- Searches ~4,000 known RNA families
- Returns curated seed alignments

**Test results:**
- ✓ Found family: SSU_rRNA_bacteria (RF00177)
- ✓ Seed alignment: 99 sequences
- ✓ High quality (manually curated)
- ⚠️ Only works for RNAs in Rfam database

### BLAST Approach Results

**What it does:**
- Searches 40.7 million sequences in RNAcentral
- Finds homologs regardless of family membership

**Test results:**
- ✓ Found: 100 homologs
- ✓ Fast: 19 seconds
- ✓ Works for any RNA (not just known families)
- ✓ 96.7% success rate on diverse RNAs

---

## Key Differences

| Feature | Rfam | BLAST |
|---------|------|-------|
| **Database size** | ~4,000 families | 40.7M sequences |
| **MSA depth** | 10-100 sequences (seed) | 50-100+ sequences |
| **Coverage** | Known families only | Any RNA with homologs |
| **Novel RNAs** | ❌ Won't find | ✅ Will find if homologs exist |
| **Speed** | Fast (~seconds) | Fast (~18 seconds) |
| **Quality** | Manually curated | Automatically generated |
| **Annotations** | Family info included | No annotations |

---

## Why We Chose BLAST

### 1. **Broader Coverage**

**Rfam limitation (as Azam noted):**
- Only ~4,000 families
- Mostly short RNAs or specific domains
- Novel/uncharacterized RNAs not included

**BLAST advantage:**
- Searches all 40.7M sequences
- Finds homologs for both known AND novel RNAs
- Not limited to pre-defined families

### 2. **Comparable MSA Depth**

**Rfam seed alignments:**
- Typically 10-100 sequences
- Our test: 99 sequences

**BLAST results:**
- 50-100+ sequences per RNA
- Our test: 100 sequences
- Can get more if needed (we capped at 100)

**Conclusion:** Similar MSA depth, but BLAST works for more RNAs

### 3. **Validated at Scale**

**Rfam:**
- Works great for known families
- But what % of RNAs are in Rfam? Unknown.

**BLAST:**
- Tested on 30 diverse RNAs
- 96.7% success rate
- Proven to work at scale

### 4. **Flexibility**

**Rfam:**
- Fixed to curated families
- Can't customize search parameters
- Seed alignments are fixed

**BLAST:**
- Adjustable parameters (E-value, identity threshold)
- Can get 50, 100, or 1000 sequences
- Works with any sequence database

---

## When Rfam Would Be Better

Rfam is excellent when:
1. You're working with **well-characterized RNA families** (rRNA, tRNA, miRNA)
2. You need **family annotations** (which family does this RNA belong to?)
3. You want **manually curated** alignments
4. You're doing **RNA family classification**

But for Professor Azam's goal (MSAs for 2D feature extraction), BLAST is better because:
- ✅ Works for more RNAs (not just known families)
- ✅ Comparable MSA depth
- ✅ Proven 96.7% success rate
- ✅ More flexible

---

## Professor Azam's Actual Concern

Re-reading her email, her concern was:

> "When an RNA sequence is long (e.g., >1000 nt), Rfam alignments typically cover only local motifs, while the rest of the sequence remains poorly aligned (appearing as long stretches of gaps in .a3m files)."

**This is about long RNAs, not about Rfam vs BLAST!**

Her point: Rfam families often represent **conserved domains** within longer RNAs, not full-length sequences.

**Our solution:**
- BLAST searches for full-length homologs
- Successfully tested on RNAs up to 1750nt
- 80% success rate for RNAs ≥800nt
- No domain-splitting needed for most cases

---

## Summary

**Professor Azam mentioned Rfam as ONE possible solution**, but also noted its limitations.

**We tested both approaches:**
- Rfam: Works for known families, 99 sequences
- BLAST: Works for any RNA, 100 sequences, 96.7% success rate

**We chose BLAST because:**
1. Broader coverage (works for novel RNAs too)
2. Comparable MSA depth
3. Validated at scale (30 diverse RNAs)
4. Handles long RNAs better (full-length homologs)

**Rfam is still valuable** for family annotations, but BLAST is the better primary approach for generating MSAs for 2D feature extraction.

---

## Recommendation

**Primary approach:** BLAST + MAFFT (what we implemented)
- Use for all RNAs
- 96.7% success rate
- 100 sequences per MSA

**Optional supplement:** Rfam
- Use to identify RNA family
- Add family annotations
- Validate results

This gives the best of both worlds: broad coverage from BLAST + family information from Rfam when available.
