import argparse
import utils
import copy

SUB = 0
CHILDREN = 1


def get_args():
    parser = argparse.ArgumentParser(description="Suffix Tree")

    parser.add_argument("--reference", help="Reference sequence file", type=str)

    parser.add_argument("--string", help="Reference sequence", type=str)

    parser.add_argument("--query", help="Query sequences", nargs="+", type=str)

    return parser.parse_args()


def add_suffix(nodes, suf):
    n = 0
    i = 0
    while i < len(suf):
        b = suf[i]
        children = nodes[n][CHILDREN]
        if b not in children:
            n2 = len(nodes)
            nodes.append([suf[i:], {}])
            nodes[n][CHILDREN][b] = n2
            return
        else:
            n2 = children[b]

        sub2 = nodes[n2][SUB]
        j = 0
        while j < len(sub2) and i + j < len(suf) and suf[i + j] == sub2[j]:
            j += 1

        if j < len(sub2):
            n3 = n2
            n2 = len(nodes)
            nodes.append([sub2[:j], {sub2[j]: n3}])
            nodes[n3][SUB] = sub2[j:]
            nodes[n][CHILDREN][b] = n2

        i += j
        n = n2


def build_suffix_tree(text):
    if text is None or text == "":
        return None
    text += "$"

    nodes = [["", {}]]

    for i in range(len(text)):
        add_suffix(nodes, text[i:])

    return nodes


def match_from_node(tree, p, node_idx, start_pos):
    if start_pos >= len(p):
        return start_pos

    node = tree[node_idx]
    children = node[CHILDREN]

    char = p[start_pos]
    if char not in children:
        return start_pos

    next_node_idx = children[char]
    next_node = tree[next_node_idx]
    substr = next_node[SUB]

    j = 0
    while j < len(substr) and start_pos + j < len(p) and substr[j] == p[start_pos + j]:
        j += 1

    if j == len(substr):
        return match_from_node(tree, p, next_node_idx, start_pos + j)
    else:
        return start_pos + j


def search_tree(suffix_tree, p):
    if p is None or p == "" or suffix_tree is None:
        return 0

    max_match = 0
    p2 = copy.deepcopy(p)

    for i in range(len(p)):
        match_len = match_from_node(suffix_tree, p, 0, 0)

        if match_len == len(p):
            max_match = match_len
            break

        p = p[1:]

    for i in range(len(p2)):
        match_len = match_from_node(suffix_tree, p2, 0, 0)

        if match_len == len(p2) and match_len > max_match:
            return match_len
        elif match_len == len(p2):
            return max_match

        p2 = p2[:-1]

    return 0


def main():
    args = get_args()

    T = None

    if args.string:
        T = args.string
    elif args.reference:
        reference = utils.read_fasta(args.reference)
        T = reference[0][1]

    tree = build_suffix_tree(T)

    if args.query:
        for query in args.query:
            match_len = search_tree(tree, query)
            print(f"{query} : {match_len}")


if __name__ == "__main__":
    main()
