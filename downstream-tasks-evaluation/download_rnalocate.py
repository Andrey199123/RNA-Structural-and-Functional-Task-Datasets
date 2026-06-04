import urllib.request
import os
from pathlib import Path

base_url = "http://www.rna-society.org/rnalocate/download/"

datasets = {
    "all_rna": "all_RNA_subcellular_localization.txt",
    "mrna": "mRNA_subcellular_localization.txt", 
    "lncrna": "lncRNA_subcellular_localization.txt",
    "circrna": "circRNA_subcellular_localization.txt",
    "mirna": "miRNA_subcellular_localization.txt",
    "pirna": "piRNA_subcellular_localization.txt",
    "trna": "tRNA_subcellular_localization.txt",
    "rrna": "rRNA_subcellular_localization.txt",
    "v2.0": "RNALocate_v2.0.txt",
}

output_dir = Path("data/raw/task2_subcellular")
output_dir.mkdir(parents=True, exist_ok=True)

print("downloading rnalocate datasets")
print("=" * 50)

for name, filename in datasets.items():
    url = base_url + filename
    output_path = output_dir / filename
    
    if output_path.exists():
        print(f"already exists: {filename}")
        continue
    
    print(f"downloading {name}: {filename}")
    
    try:
        urllib.request.urlretrieve(url, output_path)
        file_size = output_path.stat().st_size / (1024 * 1024)
        print(f"saved {filename} ({file_size:.2f} MB)")
    except Exception as e:
        print(f"failed to download {filename}: {e}")

print("\n" + "=" * 50)
print("download complete")
print(f"files saved to: {output_dir.absolute()}")
