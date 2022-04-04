from pysdd.sdd import SddManager, Vtree
import time

def identity_matrix(x,y):
    return (~x & ~y) | (x & y)

def hadamard_matrix(x,y):
    return ~y | y & ~x

def hadamard_n(start, end, x_vars, y_vars):
    if start == end - 1 :
        return hadamard_matrix(x_vars[start], y_vars[start])
    mid = int((end - start)/2 + start)
    return hadamard_n(start, mid, x_vars, y_vars) & hadamard_n(mid, end, x_vars, y_vars)


def identity_n(start, end, x_vars, y_vars):
    if start == end - 1 :
        return identity_matrix(x_vars[start], y_vars[start])
    mid = int((end - start)/2 + start)
    return identity_n(start, mid, x_vars, y_vars) & identity_n(mid, end, x_vars, y_vars)



n = pow(2, 20)
var_count = n
var_order = []
for i in range(1, n+1):
    var_order.append(i)

vtree_type = 'balanced'
vtree = Vtree(var_count, var_order, vtree_type)
manager = SddManager.from_vtree(vtree)
x_vars = []
y_vars = []
for i in range(1, n+1,2):
    x_vars.append(manager.literal(i))
    y_vars.append(manager.literal(i+1))

H = hadamard_n(0, int(n/2), x_vars, y_vars)
I = identity_n(0, int(n/2), x_vars, y_vars)

start = time.time()
ans = H + I
for i in range(1, 1024):
    ans = ans + H + I


ans.ref()
manager.minimize()
end = time.time()
print("Duration: ", end - start)
print("Size: ", ans.size())
