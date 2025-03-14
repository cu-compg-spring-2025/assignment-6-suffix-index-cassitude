import argparse
import numpy as np
import utils
import tracemalloc
import time
import random
import matplotlib.pyplot as plt
import suffix_trie
import suffix_tree
import suffix_array


def get_args():
    parser = argparse.ArgumentParser(description="Simulations")

    parser.add_argument("--reference", help="Reference sequence file", type=str)

    parser.add_argument(
        "--ref_length", help="Length of reference sequence", nargs="+", type=int
    )

    parser.add_argument(
        "--n_reads", help="Number of reads to simulate", type=int, default=5
    )

    parser.add_argument(
        "--error_rate", help="Error rate for simulated reads", type=float, default=0.05
    )

    parser.add_argument(
        "--n_size", help="Size of read lengths to simulate", type=int, default=200
    )

    return parser.parse_args()


def run_test(test_function, T, P):
    start = time.monotonic_ns()
    tracemalloc.start()
    r = test_function(T, P)
    mem = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    stop = time.monotonic_ns()

    return stop - start, mem[1] - mem[0]


def main():
    args = get_args()

    min = args.ref_length[0]
    max = args.ref_length[1]
    step = args.ref_length[2]

    trie_time_avg = []
    trie_mem_avg = []
    tree_time_avg = []
    tree_mem_avg = []
    array_time_avg = []
    array_mem_avg = []

    for r in np.arange(min, max, step):
        print(r)
        reference = utils.read_fasta(args.reference)[0][1]
        reference = reference[random.randint(0, len(reference) - r) :]
        reference = reference[:r]

        reads = utils.sim_reads(reference, args.n_size, args.n_reads, args.error_rate)

        start = time.monotonic_ns()
        tracemalloc.start()
        trie = suffix_trie.build_suffix_trie(reference)
        mem = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        stop = time.monotonic_ns()

        amor_trie_time = stop - start
        amor_trie_mem = mem[1] - mem[0]

        trie_time = []
        trie_mem = []
        for read in reads:
            total_time, mem = run_test(suffix_trie.search_trie, trie, read)
            trie_time.append(total_time)
            trie_mem.append(mem)

        start = time.monotonic_ns()
        tracemalloc.start()
        tree = suffix_tree.build_suffix_tree(reference)
        mem = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        stop = time.monotonic_ns()

        amor_tree_time = stop - start
        amor_tree_mem = mem[1] - mem[0]

        tree_time = []
        tree_mem = []
        for read in reads:
            total_time, mem = run_test(suffix_tree.search_tree, tree, read)
            tree_time.append(total_time)
            tree_mem.append(mem)

        start = time.monotonic_ns()
        tracemalloc.start()
        array = suffix_array.build_suffix_array(reference)
        mem = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        stop = time.monotonic_ns()
        amor_array_time = stop - start
        amor_array_mem = mem[1] - mem[0]

        array_time = []
        array_mem = []
        for read in reads:
            start = time.monotonic_ns()
            tracemalloc.start()
            suffix_array.search_array(reference, array, reads)
            mem = tracemalloc.get_traced_memory()
            tracemalloc.stop()
            stop = time.monotonic_ns()

            total_time = stop - start
            mem = mem[1] - mem[0]

            array_time.append(total_time)
            array_mem.append(mem)
        trie_time_avg.append(np.mean(trie_time) + amor_trie_time)
        trie_mem_avg.append(np.mean(trie_mem) + amor_trie_mem)
        tree_time_avg.append(np.mean(tree_time) + amor_tree_time)
        tree_mem_avg.append(np.mean(tree_mem) + amor_tree_mem)
        array_time_avg.append(np.mean(array_time) + amor_array_time)
        array_mem_avg.append(np.mean(array_mem) + amor_array_mem)
    # make plots
    plt.figure()
    plt.plot(
        np.arange(min, max, step),
        np.array(trie_time_avg),
        label="Trie",
    )
    plt.plot(
        np.arange(min, max, step),
        np.array(tree_time_avg),
        label="Tree",
    )
    plt.plot(np.arange(min, max, step), array_time_avg, label="Array")
    plt.xlabel("Reference length")
    plt.ylabel("Time (ns)")
    plt.legend()
    plt.savefig("figs/time.png")

    plt.figure()
    plt.plot(np.arange(min, max, step), np.array(trie_mem_avg), label="Trie")
    plt.plot(np.arange(min, max, step), np.array(tree_mem_avg), label="Tree")
    plt.plot(np.arange(min, max, step), array_mem_avg, label="Array")
    plt.xlabel("Reference length")
    plt.ylabel("Memory (bytes)")
    plt.legend()
    plt.savefig("figs/memory.png")


if __name__ == "__main__":
    main()
