#!/usr/bin/env python3
"""Parse testLuby.txt and produce a publication-quality comparison plot."""

import re
import matplotlib
matplotlib.use('Agg')  # non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ---- Parse testLuby.txt ----
data = []
with open('testLuby.txt') as f:
    text = f.read()

blocks = re.split(r'=== n=(\d+), m=(\d+) ===', text)
# blocks[0] is empty, then triples: (n, m, body)
i = 1
while i + 2 < len(blocks):
    n = int(blocks[i])
    m = int(blocks[i+1])
    body = blocks[i+2]

    yield_len = int(re.search(r'Yield length:\s*(\d+)', body).group(1))
    seq_rules = int(re.search(r'Rules:\s*(\d+)', body).group(1))
    seq_rhs = int(re.search(r'Total RHS symbols:\s*(\d+)', body).group(1))
    add_nodes = int(re.search(r'ADD nodes:\s*(\d+)', body).group(1))
    add_vars = int(re.search(r'Variables used:\s*(\d+)', body).group(1))

    cf_match = re.search(r'=== CFLOBDD ===.*?Nodes:\s*(\d+).*?Edges:\s*(\d+).*?Return map entries:\s*(\d+).*?Total:\s*(\d+)', body, re.DOTALL)
    if cf_match:
        cf_nodes = int(cf_match.group(1))
        cf_edges = int(cf_match.group(2))
        cf_retmap = int(cf_match.group(3))
        cf_total = int(cf_match.group(4))
    else:
        cf_nodes = cf_edges = cf_retmap = cf_total = 0

    data.append({
        'n': n, 'm': m, 'yield': yield_len,
        'seq_rules': seq_rules, 'seq_rhs': seq_rhs,
        'add_nodes': add_nodes, 'add_vars': add_vars,
        'cf_nodes': cf_nodes, 'cf_edges': cf_edges,
        'cf_retmap': cf_retmap, 'cf_total': cf_total,
    })
    i += 3

ns = [d['n'] for d in data]
yields = [d['yield'] for d in data]
seq_rhs = [d['seq_rhs'] for d in data]
add_nodes = [d['add_nodes'] for d in data]
cf_total = [d['cf_total'] for d in data]
cf_nodes = [d['cf_nodes'] for d in data]

# ---- Plot ----
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Left panel: all representations on log scale
ax1.semilogy(ns, yields, 'k--', linewidth=1, alpha=0.5, label='Yield length (string size)')
ax1.semilogy(ns, seq_rhs, 's-', color='#2196F3', markersize=5, linewidth=2, label='SEQUITUR RHS symbols')
ax1.semilogy(ns, add_nodes, 'D-', color='#FF9800', markersize=5, linewidth=2, label='ADD nodes')
ax1.semilogy(ns, cf_total, 'o-', color='#4CAF50', markersize=5, linewidth=2, label='CFLOBDD total (nodes+edges)')
ax1.semilogy(ns, cf_nodes, '^-', color='#9C27B0', markersize=5, linewidth=2, label='CFLOBDD nodes only')

ax1.set_xlabel('n (Luby sequence parameter)', fontsize=12)
ax1.set_ylabel('Size (log scale)', fontsize=12)
ax1.set_title('Luby Sequence: Representation Size vs n', fontsize=13)
ax1.legend(fontsize=9, loc='upper left')
ax1.grid(True, which='both', alpha=0.3)
ax1.set_xlim(-0.5, 25.5)
ax1.set_xticks(range(0, 26, 5))

# Right panel: linear scale for the three compressed representations
ax2.plot(ns, seq_rhs, 's-', color='#2196F3', markersize=5, linewidth=2, label='SEQUITUR RHS symbols')
ax2.plot(ns, add_nodes, 'D-', color='#FF9800', markersize=5, linewidth=2, label='ADD nodes')
ax2.plot(ns, cf_total, 'o-', color='#4CAF50', markersize=5, linewidth=2, label='CFLOBDD total (nodes+edges)')
ax2.plot(ns, cf_nodes, '^-', color='#9C27B0', markersize=5, linewidth=2, label='CFLOBDD nodes only')

ax2.set_xlabel('n (Luby sequence parameter)', fontsize=12)
ax2.set_ylabel('Size (linear scale)', fontsize=12)
ax2.set_title('Compressed Representations (linear)', fontsize=13)
ax2.legend(fontsize=9, loc='upper left')
ax2.grid(True, alpha=0.3)
ax2.set_xlim(-0.5, 25.5)
ax2.set_xticks(range(0, 26, 5))

plt.tight_layout()
plt.savefig('luby_comparison.png', dpi=150, bbox_inches='tight')
plt.savefig('luby_comparison.pdf', bbox_inches='tight')
print("Saved luby_comparison.png and luby_comparison.pdf")

# ---- Print table ----
print("\n{:>3s}  {:>12s}  {:>8s}  {:>9s}  {:>9s}  {:>12s}".format(
    'n', 'Yield', 'SEQ RHS', 'ADD nodes', 'CF nodes', 'CF total'))
print('-' * 60)
for d in data:
    print("{:3d}  {:12,d}  {:8d}  {:9d}  {:9d}  {:12d}".format(
        d['n'], d['yield'], d['seq_rhs'], d['add_nodes'],
        d['cf_nodes'], d['cf_total']))
