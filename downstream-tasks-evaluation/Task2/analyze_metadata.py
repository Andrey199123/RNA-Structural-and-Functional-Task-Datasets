#!/usr/bin/env python3
import argparse
import csv
import os
import json
from collections import Counter

def main():
    p = argparse.ArgumentParser(description='Analyze RNALocate metadata file')
    p.add_argument('--input', required=True)
    p.add_argument('--outdir', required=True)
    p.add_argument('--sample', type=int, default=1000)
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    counters = {
        'RNA_Type': Counter(),
        'Subcellular_Localization': Counter(),
        'Species': Counter(),
    }
    total = 0
    missing = 0
    samples = []

    with open(args.input, 'r', encoding='utf-8', errors='replace') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            total += 1
            try:
                rtype = row.get('RNA_Type','').strip()
                loc = row.get('Subcellular_Localization','').strip()
                species = row.get('Species','').strip()
            except Exception:
                missing += 1
                continue
            if not (rtype or loc or species):
                missing += 1
            counters['RNA_Type'][rtype] += 1
            counters['Subcellular_Localization'][loc] += 1
            counters['Species'][species] += 1

            if len(samples) < args.sample:
                samples.append(row)

            if total % 200000 == 0:
                print(f'Processed {total} lines...')

    summary = {
        'total_rows': total,
        'rows_with_missing_fields': missing,
    }

    # Write counters
    def write_counter(name, counter):
        outp = os.path.join(args.outdir, f'counts_by_{name}.csv')
        with open(outp, 'w', encoding='utf-8') as outfh:
            outfh.write(f'label,count\n')
            for k,v in counter.most_common():
                outfh.write(f'"{k}",{v}\n')

    write_counter('RNA_Type', counters['RNA_Type'])
    write_counter('Subcellular_Localization', counters['Subcellular_Localization'])
    write_counter('Species', counters['Species'])

    # Write summary json
    with open(os.path.join(args.outdir, 'metadata_summary.json'), 'w', encoding='utf-8') as jfh:
        json.dump({'summary': summary, 'top_localizations': counters['Subcellular_Localization'].most_common(50)}, jfh, indent=2)

    # Write sample rows
    sample_out = os.path.join(args.outdir, 'sample_rows.tsv')
    if samples:
        keys = samples[0].keys()
        with open(sample_out, 'w', encoding='utf-8') as sfh:
            writer = csv.DictWriter(sfh, fieldnames=keys, delimiter='\t')
            writer.writeheader()
            for r in samples:
                writer.writerow(r)

    print('Done. Wrote outputs to', args.outdir)

if __name__ == '__main__':
    main()
