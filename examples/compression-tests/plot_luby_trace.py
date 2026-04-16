#!/usr/bin/env python3
"""Parse testLubyTrace.txt and produce a publication-quality comparison plot."""

import re
import matplotlib
matplotlib.use('Agg')  # non-interactive backend
import matplotlib.pyplot as plt

# ---- Parse testLubyTrace.txt ----
data = []
with open('testLubyTrace.txt') as f:
    text = f.read()

blocks = re.split(r'=== n=(\d+), m=(\d+) ===', text)
i = 1
while i + 2 < len(blocks):
    n = int(blocks[i])
    m = int(blocks[i+1])
    body = blocks[i+2]

    trace_len = int(re.search(r'Trace length:\s*(\d+)', body).group(1))
    num_vars = int(re.search(r'Variables:\s*(\d+)', body).group(1))
    add_nodes = int(re.search(r'ADD nodes:\s*(\d+)', body).group(1))

    cf_match = re.search(r'=== CFLOBDD ===.*?Nodes:\s*(\d+).*?Edges:\s*(\d+).*?Return map entries:\s*(\d+).*?Total:\s*(\d+)', body, re.DOTALL)
    if cf_match:
        cf_nodes = int(cf_match.group(1))
        cf_edges = int(cf_match.group(2))
        cf_retmap = int(cf_match.group(3))
        cf_total = int(cf_match.group(4))
    else:
        cf_nodes = cf_edges = cf_retmap = cf_total = 0

    data.append({
        'n': n, 'm': m, 'trace_len': trace_len,
        'num_vars': num_vars, 'add_nodes': add_nodes,
        'cf_nodes': cf_nodes, 'cf_edges': cf_edges,
        'cf_retmap': cf_retmap, 'cf_total': cf_total,
    })
    i += 3

# Also load SEQUITUR data for comparison (if available)
seq_data = []
try:
    with open('testLuby.txt') as f:
        text = f.read()
    blocks = re.split(r'=== n=(\d+), m=(\d+) ===', text)
    i = 1
    while i + 2 < len(blocks):
        n = int(blocks[i])
        body = blocks[i+2]
        seq_rhs = int(re.search(r'Total RHS symbols:\s*(\d+)', body).group(1))
        add_nodes = int(re.search(r'ADD nodes:\s*(\d+)', body).group(1))
        cf_match = re.search(r'=== CFLOBDD ===.*?Total:\s*(\d+)', body, re.DOTALL)
        cf_total = int(cf_match.group(1)) if cf_match else 0
        seq_data.append({'n': n, 'seq_rhs': seq_rhs, 'seq_add': add_nodes, 'seq_cf_total': cf_total})
        i += 3
except FileNotFoundError:
    pass

ns = [d['n'] for d in data]
trace_lens = [d['trace_len'] for d in data]
add_nodes = [d['add_nodes'] for d in data]
num_vars = [d['num_vars'] for d in data]
cf_total = [d['cf_total'] for d in data]
cf_nodes = [d['cf_nodes'] for d in data]

# ---- Plot ----
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle('Luby Sequence: Direct Trace→ADD vs SEQUITUR→ADD', fontsize=14)

# Left panel: log scale
ax1.semilogy(ns, trace_lens, 'k--', linewidth=1, alpha=0.5, label='Trace length')
ax1.semilogy(ns, add_nodes, 'D-', color='#FF9800', markersize=5, linewidth=2, label='Direct ADD nodes')
ax1.semilogy(ns, cf_total, 'o-', color='#4CAF50', markersize=5, linewidth=2, label='Direct CFLOBDD total')
ax1.semilogy(ns, cf_nodes, '^-', color='#9C27B0', markersize=5, linewidth=2, label='Direct CFLOBDD nodes')
if seq_data:
    sns = [d['n'] for d in seq_data]
    ax1.semilogy(sns, [d['seq_rhs'] for d in seq_data], 's--', color='#2196F3', markersize=4, linewidth=1.5, alpha=0.7, label='SEQUITUR RHS')
    ax1.semilogy(sns, [d['seq_add'] for d in seq_data], 'D--', color='#FF9800', markersize=4, linewidth=1.5, alpha=0.4, label='SEQUITUR ADD nodes')
    ax1.semilogy(sns, [d['seq_cf_total'] for d in seq_data], 'o--', color='#4CAF50', markersize=4, linewidth=1.5, alpha=0.4, label='SEQUITUR CFLOBDD total')

ax1.set_xlabel('n (Luby sequence parameter)', fontsize=12)
ax1.set_ylabel('Size (log scale)', fontsize=12)
ax1.set_title('All representations (log scale)', fontsize=12)
ax1.legend(fontsize=8, loc='upper left')
ax1.grid(True, which='both', alpha=0.3)
ax1.set_xlim(-0.5, 25.5)
ax1.set_xticks(range(0, 26, 5))

# Right panel: linear scale, compressed representations
ax2.plot(ns, add_nodes, 'D-', color='#FF9800', markersize=5, linewidth=2, label='Direct ADD nodes')
ax2.plot(ns, cf_total, 'o-', color='#4CAF50', markersize=5, linewidth=2, label='Direct CFLOBDD total')
ax2.plot(ns, cf_nodes, '^-', color='#9C27B0', markersize=5, linewidth=2, label='Direct CFLOBDD nodes')
ax2.plot(ns, num_vars, 'x-', color='#607D8B', markersize=5, linewidth=1.5, label='Direct ADD variables')
if seq_data:
    ax2.plot(sns, [d['seq_rhs'] for d in seq_data], 's--', color='#2196F3', markersize=4, linewidth=1.5, alpha=0.7, label='SEQUITUR RHS')
    ax2.plot(sns, [d['seq_add'] for d in seq_data], 'D--', color='#FF9800', markersize=4, linewidth=1.5, alpha=0.4, label='SEQUITUR ADD nodes')
    ax2.plot(sns, [d['seq_cf_total'] for d in seq_data], 'o--', color='#4CAF50', markersize=4, linewidth=1.5, alpha=0.4, label='SEQUITUR CFLOBDD total')

ax2.set_xlabel('n (Luby sequence parameter)', fontsize=12)
ax2.set_ylabel('Size (linear scale)', fontsize=12)
ax2.set_title('Compressed representations (linear)', fontsize=12)
ax2.legend(fontsize=8, loc='upper left')
ax2.grid(True, alpha=0.3)
ax2.set_xlim(-0.5, 25.5)
ax2.set_xticks(range(0, 26, 5))

plt.tight_layout()
plt.savefig('luby_trace_comparison.png', dpi=150, bbox_inches='tight')
plt.savefig('luby_trace_comparison.pdf', bbox_inches='tight')
print("Saved luby_trace_comparison.png and luby_trace_comparison.pdf")

# ---- Print table ----
print("\n{:>3s}  {:>12s}  {:>5s}  {:>9s}  {:>9s}  {:>9s}".format(
    'n', 'Trace len', 'Vars', 'ADD nodes', 'CF nodes', 'CF total'))
print('-' * 55)
for d in data:
    print("{:3d}  {:12,d}  {:5d}  {:9d}  {:9d}  {:9d}".format(
        d['n'], d['trace_len'], d['num_vars'], d['add_nodes'],
        d['cf_nodes'], d['cf_total']))
