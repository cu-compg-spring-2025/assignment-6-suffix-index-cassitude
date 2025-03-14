import argparse
import utils


def get_args():
    parser = argparse.ArgumentParser(description="Suffix Trie")

    parser.add_argument("--reference", help="Reference sequence file", type=str)

    parser.add_argument("--string", help="Reference sequence", type=str)

    parser.add_argument("--query", help="Query sequences", nargs="+", type=str)

    return parser.parse_args()


def build_suffix_trie(s):
    if s is None or s == "":
        return None

    trie = {}
    suffs = [s[i:] for i in range(len(s))]

    for suffix in suffs:
        current_node = trie
        for char in suffix:
            if char not in current_node:
                current_node[char] = {}
            current_node = current_node[char]
        current_node["$"] = {}

    return trie


def search_trie(trie, p):
    if trie is None or p is None or p == "":
        return 0

    current_node = trie
    match_length = 0

    for char in p:
        if char in current_node:
            match_length += 1
            current_node = current_node[char]
        else:
            return match_length

    return match_length


def main():
    args = get_args()

    T = None

    if args.string:
        T = args.string
    elif args.reference:
        reference = utils.read_fasta(args.reference)
        T = reference[0][1]

    trie = build_suffix_trie(T)

    if args.query:
        for query in args.query:
            match_len = search_trie(trie, query)
            print(f"{query} : {match_len}")


if __name__ == "__main__":
    main()
