import os
import sys
import subprocess
from pathlib import Path

def run_command(cmd, cwd=None):
    try:
        result = subprocess.run(cmd, shell=True, cwd=cwd, capture_output=True, text=True)
        if result.returncode == 0:
            return True, result.stdout
        else:
            return False, result.stderr
    except Exception as e:
        return False, str(e)



def clone_git_repo(repo_url, output_dir):
    output_dir = Path(output_dir)
    
    if output_dir.exists():
        print(f"repository already exists: {output_dir}")
        return True
    
    print(f"cloning {repo_url}")
    
    success, output = run_command(f"git clone {repo_url} {output_dir}")
    if success:
        print(f"cloned to {output_dir}")
        return True
    else:
        print(f"clone failed: {output}")
        return False

def download_task1_rsa():
    print("\n=== task 1: rna solvent accessibility (rsa) ===")
    
    base_dir = Path("data/raw/task1_rsa")
    base_dir.mkdir(parents=True, exist_ok=True)
    
    repo_url = "https://github.com/jaswindersingh2/RNAsnap2.git"
    repo_dir = base_dir / "RNAsnap2"
    
    return clone_git_repo(repo_url, repo_dir)

def download_task2_subcellular():
    print("\n=== task 2: rna subcellular localization (under investigation) ===")
    
    base_dir = Path("data/raw/task2_subcellular")
    base_dir.mkdir(parents=True, exist_ok=True)
    
    print("note: rnalocate v2.0 requires manual download from rnalocate.org")
    print("visit: http://www.rna-society.org/rnalocate/")
    print("download the dataset and place it in:", base_dir.absolute())
    
    return True

def download_task3_contact_map():
    print("\n=== task 3: rna contact map prediction ===")
    
    base_dir = Path("data/raw/task3_contact_map")
    base_dir.mkdir(parents=True, exist_ok=True)
    
    print("note: rnacontact dataset requires accessing the paper supplementary materials")
    print("paper: sun et al., bioinformatics 2021")
    print("check github or paper supplementary for dataset access")
    print("expected: 301 filtered sequences from 1,786 pdb sequences")
    
    return True

def download_task4_rbp():
    print("\n=== task 4: rna-protein binding prediction ===")
    
    base_dir = Path("data/raw/task4_rbp")
    base_dir.mkdir(parents=True, exist_ok=True)
    
    repo_url = "https://github.com/kuixu/PrismNet.git"
    repo_dir = base_dir / "PrismNet"
    
    return clone_git_repo(repo_url, repo_dir)

def download_task5_rri():
    print("\n=== task 5: rna-rna interaction ===")
    
    base_dir = Path("data/raw/task5_rri")
    base_dir.mkdir(parents=True, exist_ok=True)
    
    print("note: deepmirtar dataset requires accessing the paper or contacting authors")
    print("paper: wen et al., bioinformatics 2018")
    print("expected: 27,720 mirna-mrna pairs (13,860 positive + 13,860 negative)")
    
    return True

def download_task6_secondary_structure():
    print("\n=== task 6: rna secondary structure prediction (on hold) ===")
    
    print("status: on hold - awaiting approval")
    print("do not proceed with download until explicit approval")
    
    return True

def download_task7_long_range():
    print("\n=== task 7: long-range base pair prediction ===")
    
    base_dir = Path("data/raw/task7_long_range")
    base_dir.mkdir(parents=True, exist_ok=True)
    
    print("dataset a: pdb benchmark (532 sequences)")
    print("dataset b: bprna-1m tr0 (10,814 sequences)")
    print("note: requires accessing rna-km paper supplementary materials")
    
    return True

def main():
    print("downstream tasks dataset downloader")
    print("=" * 50)
    
    os.makedirs("data/raw", exist_ok=True)
    
    tasks = [
        ("task 1: rsa", download_task1_rsa),
        ("task 2: subcellular", download_task2_subcellular),
        ("task 3: contact map", download_task3_contact_map),
        ("task 4: rbp", download_task4_rbp),
        ("task 5: rri", download_task5_rri),
        ("task 6: secondary structure", download_task6_secondary_structure),
        ("task 7: long-range", download_task7_long_range),
    ]
    
    results = {}
    
    for task_name, task_func in tasks:
        try:
            success = task_func()
            results[task_name] = "success" if success else "failed"
        except Exception as e:
            print(f"error in {task_name}: {e}")
            results[task_name] = "error"
    
    print("\n" + "=" * 50)
    print("download summary:")
    for task_name, status in results.items():
        print(f"  {task_name}: {status}")

if __name__ == "__main__":
    main()
