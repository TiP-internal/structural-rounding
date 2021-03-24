
#include "sr_apx/domset/exact/domset_exact.hpp"
#include "sr_apx/util/util.hpp"

namespace sr_apx {
namespace domset {
namespace exact {

bool intersect(const Set& a, const Set& b) {
    for (int x : a) {
        if (b.contains(x)) {
            return true;
        }
    }

    return false;
}

struct table_entry {
    int size_below = 0;
    int left_ref = -1;
    int right_ref = -1;
};

using table = Map<table_entry>;

// key format:
// k bits, 1 shows in domset
// up to k bits, 1 shows dominated

// join format:
// k bits, select
// k bits, 1-1 shows dominated left
//         1-0 shows dominated right
//         0-1 shows in domset
//         0-0 shows may be undominated

// TODO: (join) make sure tree decomp can't have equal bag as left child
// TODO: (intro) make sure always a forget from left child

// loop over parent entries, convert to child
// width wrt parent
table forget(const table& child, int index, int width, int default_size, bool copy_ref, bool opt) {
    table parent;

    int mask = (1 << index) - 1;

    util::GrayCode start;
    util::GrayCode end;
    std::tie(start, end) = util::graycode(width);

    std::vector<int> bitsinlower;
    for (int i = 0; i < width; ++i) {
        bitsinlower.push_back(i);
    }

    while (start != end) {
        util::GrayCode ls;
        util::GrayCode le;
        std::tie(ls, le) = util::graycode(bitsinlower.size());

        int lower = 0;

        while (ls != le) {
            table_entry parent_entry;
            parent_entry.size_below = default_size;

            int child_upper = ((start.get_code() & ~mask) << 1) + (start.get_code() & mask);
            int child_lower = ((lower & ~mask) << 1) + (lower & mask);

            int temp = child_upper + (1 << index);
            const table_entry& in_ds = child.at((temp << 16) + child_lower);
            if (in_ds.size_below + 1 < parent_entry.size_below) {
                parent_entry.size_below = in_ds.size_below + 1;
                if (copy_ref) {
                    parent_entry.left_ref = in_ds.left_ref;
                }
                else {
                    parent_entry.left_ref = (temp << 16) + child_lower;
                }
            }

            temp = child_lower + (1 << index);
            if (opt) {
                temp = child_lower;
            }
            const table_entry& covered = child.at((child_upper << 16) + temp);
            if (covered.size_below < parent_entry.size_below) {
                parent_entry.size_below = covered.size_below;
                if (copy_ref) {
                    parent_entry.left_ref = covered.left_ref;
                }
                else {
                    parent_entry.left_ref = (child_upper << 16) + temp;
                }
            }

            parent[(start.get_code() << 16) + lower] = parent_entry;

            ++ls;
            if (ls == le) {
                break;
            }

            int change = ls.get_change();
            lower ^= 1 << (bitsinlower[change]);
        }

        ++start;
        if (start == end) {
            break;
        }

        int change = start.get_change();
        bool flag = false;
        std::vector<int>::iterator b = bitsinlower.begin();
        while (b != bitsinlower.end()) {
            if (*b == change) {
                bitsinlower.erase(b);
                flag = true;
                break;
            }

            ++b;
        }

        if (!flag) {
            bitsinlower.push_back(change);
        }
    }

    return parent;
}

// loop over domset upper child (gc)
// loop over dominated below child
// insert extra bits to get parent index
// track number of 0 below it in upper
// bag wrt to child
table introduce(const table& child, const std::vector<int>& bag, const Set& neighbors, int index, int default_size) {
    table parent;

    int mask = (1 << index) - 1;

    util::GrayCode start;
    util::GrayCode end;
    std::tie(start, end) = util::graycode(bag.size());

    Set doms;
    std::vector<int> bitsinlower;
    for (int i = 0; i < bag.size(); ++i) {
        bitsinlower.push_back(i);
    }

    while (start != end) {
        util::GrayCode ls;
        util::GrayCode le;
        std::tie(ls, le) = util::graycode(bitsinlower.size());

        int lower = 0;
        int lower_alt = 0;

        while (ls != le) {
            const table_entry& child_entry = child.at((start.get_code() << 16) + lower);

            int parent_upper = ((start.get_code() & ~mask) << 1) + (start.get_code() & mask);
            int parent_lower = ((lower & ~mask) << 1) + (lower & mask);

            const table_entry& child_entry_alt = child.at((start.get_code() << 16) + lower_alt);
            int temp = parent_upper + (1 << index);
            table_entry in_ds;
            in_ds.size_below = child_entry_alt.size_below;
            in_ds.left_ref = child_entry_alt.left_ref;
            parent[(temp << 16) + parent_lower] = in_ds;

            temp = parent_lower + (1 << index);
            table_entry covered;
            covered.size_below = default_size;
            if (intersect(doms, neighbors)) {
                covered.size_below = child_entry.size_below;
                covered.left_ref = child_entry.left_ref;
            }
            parent[(parent_upper << 16) + temp] = covered;

            table_entry maybe;
            maybe.size_below = child_entry.size_below;
            maybe.left_ref = child_entry.left_ref;
            parent[(parent_upper << 16) + parent_lower] = maybe;

            ++ls;
            if (ls == le) {
                break;
            }

            int change = ls.get_change();
            lower ^= 1 << (bitsinlower[change]);
            if (!neighbors.contains(bag[bitsinlower[change]])) {
                lower_alt ^= 1 << (bitsinlower[change]);
            }
        }

        ++start;
        if (start == end) {
            break;
        }

        int change = start.get_change();
        bool flag = false;
        std::vector<int>::iterator b = bitsinlower.begin();
        while (b != bitsinlower.end()) {
            if (*b == change) {
                bitsinlower.erase(b);
                flag = true;
                doms.insert(bag[change]);
                break;
            }

            ++b;
        }

        if (!flag) {
            bitsinlower.push_back(change);
            doms.erase(bag[change]);
        }
    }

    return parent;
}

// loop over select
// loop over k bits
// compute left index, right index, parent index
// update size, left ref, right ref if smaller
// done without graycodes to improve memory access pattern
table join(const table& left, const table& right, int width, int default_size) {
    table parent;
    for (int upper = 0; upper < (1 << width); ++upper) {
        for (int lower = 0; lower < (1 << width); ++lower) {
            int parent_index = (((~upper) & lower) << 16) + upper;
            if (!parent.contains(parent_index)) {
                table_entry entry;
                entry.size_below = default_size;
                parent[parent_index] = entry;
            }

            table_entry& parent_entry = parent[parent_index];

            int left_index = (((~upper) & lower) << 16) + (upper & lower);
            const table_entry& left_entry = left.at(left_index);

            int right_index = (((~upper) & lower) << 16) + (upper & (~lower));
            const table_entry& right_entry = right.at(right_index);

            if (left_entry.size_below + right_entry.size_below < parent_entry.size_below) {
                parent_entry.size_below = left_entry.size_below + right_entry.size_below;
                parent_entry.left_ref = left_entry.left_ref;
                parent_entry.right_ref = right_index;
            }
        }
    }

    return parent;
}

Set tw_exact(const Graph& graph, treewidth::Decomposition& decomp, const Set& opt) {
    // ensure vertices in bags are sorted
    decomp.sort_bags();

    std::vector<table> tables(decomp.pre_order.size());

    for (int i = tables.size() - 1; i >= 0; --i) {
        std::vector<int>& parent = decomp.pre_order[i].bag;

        int child_index = decomp.pre_order[i].left_child;
        std::vector<int> default_bag;
        std::vector<int>& child = default_bag;

        table default_table;
        default_table[0] = table_entry();
        table& child_table = default_table;

        if (child_index != -1) {
            child = decomp.pre_order[child_index].bag;
            child_table = tables[child_index];
        }

        bool copy_ref = false;

        std::vector<int> use_bag;

        // forget
        int p = 0;
        int c = 0;
        int index = 0;
        int width = child.size();
        while (c < child.size()) {
            if (p == parent.size() || child[c] < parent[p]) {
                // forget child[c]
                child_table = forget(child_table, index, width - 1, graph.size() + 1, copy_ref, opt.contains(child[c]));
                copy_ref = true;
                --width;
                ++c;
            }
            else if (child[c] == parent[p]) {
                use_bag.push_back(child[c]);
                ++p;
                ++c;
                ++index;
            }
            else {
                ++p;
            }
        }

        // introduce
        p = 0;
        c = 0;
        index = 0;
        while (p < parent.size()) {
            if (c == child.size() || parent[p] < child[c]) {
                // introduce parent[p]
                child_table = introduce(child_table, use_bag, graph.neighbors(parent[p]), p, graph.size() + 1);
                use_bag.insert(use_bag.begin() + p, parent[p]);
                ++p;
            }
            else if (child[c] == parent[p]) {
                ++p;
                ++c;
            }
            else {
                ++c;
            }
        }

        child_index = decomp.pre_order[i].right_child;

        if (child_index != -1) {
            // join
            child_table = join(child_table, tables[child_index], parent.size(), graph.size() + 1);
        }

        tables[i] = child_table;
    }

    Set res;
    return res;
}

}}}
