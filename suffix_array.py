import argparse
import utils
import suffix_tree

SUB = 0
CHILDREN = 1


def get_args():
    parser = argparse.ArgumentParser(description="Suffix Tree")

    parser.add_argument("--reference", help="Reference sequence file", type=str)

    parser.add_argument("--string", help="Reference sequence", type=str)

    parser.add_argument("--query", help="Query sequences", nargs="+", type=str)

    return parser.parse_args()


def build_suffix_array(T):
    tree = suffix_tree.build_suffix_tree(T)

    suffix_array = []
    node_array = []
    visited = set()

    stack = [(0, 0)]  # (node_idx, current_path_length)
    while stack:
        node_idx, path_len = stack.pop()

        if node_idx in visited:
            continue
        visited.add(node_idx)

        node = tree[node_idx]
        substr = node[SUB]
        current_path_len = path_len + len(substr)

        if "$" in substr:
            suffix_start = current_path_len - 1
            suffix_array.append(suffix_start)
            node_array.append(node[0])

        for char, child_idx in sorted(tree[node_idx][CHILDREN].items(), reverse=True):
            stack.append((child_idx, current_path_len))
    return sorted(suffix_array, key=lambda pos: T[pos:])


def search_array(T, suffix_array, q):
    if not suffix_array or not q:
        return 0

    best_match_len = 0

    for i in range(len(q)):
        lo = 0
        hi = len(suffix_array) - 1

        while lo <= hi:
            mid = (lo + hi) // 2

            suffix_start = suffix_array[mid]
            suffix = T[suffix_start:]

            match_len = 0
            for i in range(min(len(q), len(suffix))):
                if q[i] == suffix[i]:
                    match_len += 1
                else:
                    break

            best_match_len = max(best_match_len, match_len)

            if match_len == len(q):
                return match_len

            if mid > 0 and (
                match_len < len(q)
                and (match_len == len(suffix) or q[match_len] < suffix[match_len])
            ):
                hi = mid - 1
            else:
                lo = mid + 1
        q = q[1:]

    return best_match_len


def build_suffix_array_o(T):
    tree = suffix_tree.build_suffix_tree(T)
    # Your code here

    # defs
    stack = [0]
    while stack:
        node_idx = stack.pop()
        for child in tree[node_idx][CHILDREN]:
            stack.append(child)

    return None


def search_array_o(T, suffix_array, q):
    # Your code here

    # binary search
    lo = -1
    hi = len(suffix_array)
    while hi - lo > 1:
        mid = int((lo + hi) / 2)
        if suffix_array[mid] < q:
            lo = mid
        else:
            hi = mid
    return hi


def main():
    args = get_args()

    T = None

    if args.string:
        T = args.string
    elif args.reference:
        reference = utils.read_fasta(args.reference)
        T = reference[0][1]

    array = build_suffix_array(T)

    print(array)

    if args.query:
        for query in args.query:
            match_len = search_array(T, array, query)
            print(f"{query} : {match_len}")


if __name__ == "__main__":
    main()
