#!/usr/bin/env python3
import argparse
import csv
import os
import hashlib

def assign_split(id_str, seed, train_frac, val_frac):
    key = f"{id_str}|{seed}".encode('utf-8')
    h = hashlib.md5(key).hexdigest()
    v = int(h[:8],16) / 0xFFFFFFFF
    if v < train_frac:
        return 'train'
    elif v < train_frac + val_frac:
        return 'val'
    else:
        return 'test'


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--input', required=True)
    p.add_argument('--outdir', required=True)
    p.add_argument('--seed', type=int, default=42)
    p.add_argument('--train_frac', type=float, default=0.8)
    p.add_argument('--val_frac', type=float, default=0.1)
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    inf = args.input
    out_files = {
        'train': open(os.path.join(args.outdir,'train.tsv'),'w',encoding='utf-8'),
        'val': open(os.path.join(args.outdir,'val.tsv'),'w',encoding='utf-8'),
        'test': open(os.path.join(args.outdir,'test.tsv'),'w',encoding='utf-8'),
    }
    writers = {}
    counts = {'train':0,'val':0,'test':0}

    with open(inf,'r',encoding='utf-8',errors='replace') as fh:
        reader = csv.reader(fh, delimiter='\t')
        header = next(reader)
        for k in out_files:
            out_files[k].write('\t'.join(header) + '\n')
            writers[k] = csv.writer(out_files[k], delimiter='\t')

        idx_id = None
        # find index of RNALocate_ID if present
        for i,h in enumerate(header):
            if h.strip() == 'RNALocate_ID':
                idx_id = i
                break
        if idx_id is None:
            raise SystemExit('RNALocate_ID column not found in header')

        for row in reader:
            rid = row[idx_id]
            split = assign_split(rid, args.seed, args.train_frac, args.val_frac)
            writers[split].writerow(row)
            counts[split] += 1

    for f in out_files.values():
        f.close()

    print('Split counts:', counts)

if __name__ == '__main__':
    main()
