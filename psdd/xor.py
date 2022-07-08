from pysdd.sdd import SddManager, Vtree
import time

n = pow(2, 16)
var_count = n
var_order = []
for i in range(1, n+1):
    var_order.append(i)

vtree_type = 'balanced'
vtree = Vtree(var_count, var_order, vtree_type)
manager = SddManager.from_vtree(vtree)
x_vars = []
for i in range(1, n+1):
    x_vars.append(manager.literal(i))

def xor(a,b):
    return (a & ~b) | (~a & b)

start = time.time()
f = x_vars[0]

for i in range(1, n):
    f = xor(x_vars[i], f)

f.ref()
manager.minimize()
end = time.time()
print(end - start)
print(f.size())
