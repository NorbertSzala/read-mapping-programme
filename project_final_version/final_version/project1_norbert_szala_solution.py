#!/usr/bin/env python3
# -*- coding: utf-8 -*-


#########################################
######   Import needed packages   #######
#########################################
import os
import argparse
import logging
import psutil
from Bio import SeqIO
import threading
import time

# ---- import mine functions ----
import utils.utils as ut


#########################################
########## Monitor RAM usage  ###########
#########################################


def monitor_memory(interval=0.05):
    proc = psutil.Process()
    peak = 0
    while True:
        mem = proc.memory_info().rss / (1024**2)
        peak = max(peak, mem)
        time.sleep(interval)
        if getattr(monitor_memory, "stop", False):
            break
    monitor_memory.peak = peak
    
def start_memory_monitor():
    t = threading.Thread(target=monitor_memory, daemon=True)
    t.start()
    return t

def stop_memory_monitor():
    monitor_memory.stop = True
    time.sleep(0.1)
    return monitor_memory.peak


#########################################
############# Run pipeline #############
#########################################
def run_pipeline(
    K: int = 13,                        # K-mer length. Smaller K = higher sensitivity, larger K = fewer false seeds.
    W: int = 40,                        # Window length for minimizers. Controls minimizer density (1 per W).
    cluster_gap: int = 30,             # Max allowed distance between seeds to be merged into the same cluster.
    min_cluster_seeds: int = 4,         # Minimal number of seeds for a cluster to be considered a valid candidate.
    max_cand: int = 4,                  # Max number of top-scoring clusters passed to DP verification.
    max_error_ratio: float = 0.09,      # Expected upper bound on read error rate; controls DP band width.
    reference_path: str = "",           # Path to reference FASTA file.
    reads_path: str = "",               # Path to reads FASTA file.
    sketch_size: int = 200,             # Size of bottom-k sketch used for Jaccard similarity estimation.
    ref_window_size: int = 2000,        # Size of fixed windows used for reference sketching.
    ref_stride: int = 500,              # Step between consecutive reference windows (controls resolution).
    jaccard_thresh: float | None = None # Jaccard threshold for window selection; auto-set if None.
):

    # load reference (single sequence assumed)
    ref_seq = None
    for rid, seq in ut.read_fasta_gen(reference_path):
        ref_seq = seq
        break

    if ref_seq is None:
        raise FileNotFoundError(
            f"Reference FASTA empty or not found at {reference_path}"
        )
        
    # build minimizer index 
    minimizer_index = ut.compute_minimizer_index(ref_seq, K=K, W=W)

    # precompute reference sketches for windows
    ref_windows = ut.build_reference_window_sketches(
        ref_seq,
        K=K,
        window_size=ref_window_size,
        sketch_size=sketch_size,
        stride=ref_stride,
    )
    
    # compute default jaccard threshold if not provided
    if jaccard_thresh is None:
        jaccard_thresh = max(0.02, 1.0 - max_error_ratio * 1.5)

    results = []


    # process reads
    for read_id, read in ut.read_fasta_gen(reads_path):

        # Read sketch
        read_hashes = ut.build_kmer_hashes(read, K)
        read_sketch = ut.minhash_sketch_from_kmer_hashes(read_hashes, sketch_size)

        # Compare read sketch with all reference windows
        cand_windows = []
        for (start, end, rsk) in ref_windows:
            j = ut.jaccard_from_sketches(read_sketch, rsk)
            if j >= jaccard_thresh:
                cand_windows.append((j, start, end))

        # If nothing passed threshold → fallback: take top-N windows
        if not cand_windows:
            scored = []
            for (start, end, rsk) in ref_windows:
                j = ut.jaccard_from_sketches(read_sketch, rsk)
                scored.append((j, start, end))

            scored.sort(reverse=True)
            cand_windows = scored[: max_cand]

        # Extract window coordinates
        cand_windows = [(s, e) for (_, s, e) in cand_windows]

        # --- SEEDING ---
        all_seed_hits = []
        seeds, unmapped = ut.seeding_from_minimizers(read, minimizer_index, K=K, W=W)
        if seeds:
            for (wstart, wend) in cand_windows:
                seeds_in_win = [(r, p) for (r, p) in seeds if wstart <= r < wend]
                all_seed_hits.extend(seeds_in_win)

        # --- CHAINING ---
        clusters = ut.chain_seeds(
            all_seed_hits,
            cluster_gap=cluster_gap,
            min_cluster_seeds=min_cluster_seeds,
        )

        # --- LIMIT + VERIFY ---
        cand_limited = ut.limit_candidates(clusters, max_cand=max_cand)
        verified = []
        for cand in cand_limited:
            v = ut.verify_candidate(read, ref_seq, cand, max_error_ratio=max_error_ratio)
            if v["start"] is None:
                continue
            verified.append(v)



        # --- RANKING ---
        best = ut.rank_and_output(read_id, read, verified)

        if best is None:
            results.append((read_id, "UNMAPPED"))
        else:
            results.append((read_id, best))

    return results


#########################################
################  Main   ################
#########################################

def main():

    parser = argparse.ArgumentParser(description="Read mapper pipeline")
    parser.add_argument("reference", help="Path to reference FASTA file")
    parser.add_argument("reads", help="Path to reads FASTA file")
    parser.add_argument("output", help="Path to output file")
    args = parser.parse_args()
    
    # Setup logging
    log_path = os.path.join(os.path.dirname(args.output), "pipeline.log")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        filename=log_path,
        filemode="w",
    )

    # Start total timer and memory
    start_total_time = time.perf_counter()
    mem_thread = start_memory_monitor()

    
    # Execute
    logging.info("Starting pipeline")
    results = run_pipeline(
        K=13,
        W=40,
        cluster_gap=30,
        min_cluster_seeds=4,
        max_cand=4,
        max_error_ratio=0.09,
        reference_path=args.reference,
        reads_path=args.reads,
    )


    # End total timer and memory
    end_total_time = time.perf_counter()
    num_reads = sum(1 for _ in SeqIO.parse(args.reads, "fasta"))
    total_time = end_total_time - start_total_time
    peak_mem = stop_memory_monitor()
    
    # Write output
    with open(args.output, "w") as f:
        for rid, out in results:
            if out == "UNMAPPED":
                f.write(f"{rid}\tUNMAPPED\n")
            else:
                f.write(f"{rid}\t{out[0]}\t{out[1]}\n")
                
    mapped = sum(1 for _, out in results if out != "UNMAPPED")
    total = len(results)
    rate = mapped / total if total > 0 else 0.0
    
    
    logging.info(f"Pipeline finished.")
    logging.info(f"Mapping summary: {mapped}/{total} reads mapped ({rate:.2%}).")
    logging.info(f"Total runtime: {total_time:.2f} seconds")
    logging.info(f"Approximate memory used: {peak_mem:.2f} MB")

    print(f"Pipeline finished in {total_time:.2f} seconds")
    print(f"[INFO] Mapping summary: {mapped}/{total} reads mapped ({rate:.2%}).")
    print(f'Time of processing one read: {(num_reads/total_time):.2f}')
    print(f"Approximate memory used: {peak_mem:.2f} MB")
    print(f"Results written to {args.output}")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.error("Pipeline crashed", exc_info=True)
        raise
    
