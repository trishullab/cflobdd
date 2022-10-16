import math
import time
import python_wrapper_cflobdd as pwc

def get_idx(i,j,dim, dim2=None):
    if dim2 is None:
        dim2 = dim
    if (i < dim2//2 and j < dim//2):
        return dim2//2 * i + j + 1
    elif (i < dim2//2 and j >= dim//2):
        return (dim//2 * dim2//2) + dim2//2 * i + (j - dim2//2) + 1
    elif (i >= dim2//2 and j < dim//2):
        return dim//2 * dim2 + (i - dim//2)*dim2//2 + j + 1
    else:
        return (dim//2 * dim2) + (dim//2 * dim2//2) + dim2//2 * (i - dim//2) + (j - dim2//2) + 1
    return dim2*i+j+1



def get_neighbors(i,j,dim, dim2=None):
    if dim2 is None:
        dim2 = dim
    ret = []
    d = [(-1,-1),(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0),(1,1)]
    for x,y in d:
        ii = i+x
        jj = j+y
        if ii >= 0 and jj >= 0 and ii < dim and jj < dim2:
            ret.append((ii,jj))
    return ret


def compute_alpha(betas, start=0, end = 0):
    print("beg: ", len(betas), "start: ", start, " end: ", end)
    if len(betas) == 1:
        return betas[0]
    if len(betas) == 2:
        return pwc.CFLOBDD.MkAnd(betas[0], betas[1])

    mid = len(betas)//2
    a1 = compute_alpha(betas[:mid], start, start + mid)
    a2 = compute_alpha(betas[mid:], start + mid, end)

    a3 = pwc.CFLOBDD.MkAnd(a1, a2)
    print(len(betas), "start: ", start, " end: ", end)
    return a3


# assumes output is dim*dim square
def build_constraints(dim, save_path=None):
    var_count = dim * dim
    # final constraint variables
    level = math.ceil(math.log(var_count, 2))
    alpha1 = pwc.CFLOBDD.MkTrue(level)
    alpha2 = pwc.CFLOBDD.MkTrue(level)
    alpha3 = pwc.CFLOBDD.MkTrue(level)
    alpha4 = pwc.CFLOBDD.MkTrue(level)
    var_cflobdds = [pwc.CFLOBDD.MkProjection(i, level) for i in range(var_count)]
    # list_of_betas = []

    for i in range(dim//2):
        for j in range(dim//2):
            if i == 0 and j == 0:
                continue
            else:
                nbrs = get_neighbors(i, j, dim)
                beta = pwc.CFLOBDD.MkFalse(level)
                for a in range(len(nbrs)):
                    a_idx = get_idx(nbrs[a][0], nbrs[a][1], dim)
                    gamma = pwc.CFLOBDD.MkTrue(level)
                    for b in range(len(nbrs)):
                        if a != b:
                            id_b = get_idx(nbrs[b][0], nbrs[b][1], dim)
                            gamma = pwc.CFLOBDD.MkAnd(gamma, pwc.CFLOBDD.MkNot(var_cflobdds[id_b-1]))

                    gamma = pwc.CFLOBDD.MkAnd(gamma, var_cflobdds[a_idx-1])
                    beta = pwc.CFLOBDD.MkOr(beta, gamma)

                idx = get_idx(i, j, dim)-1
                #beta = pwc.CFLOBDD.MkAnd(beta, var_cflobdds[idx])
                beta = pwc.CFLOBDD.MkOr(beta, pwc.CFLOBDD.MkNot(var_cflobdds[idx]))
                #list_of_betas.append(beta)
                alpha1 = pwc.CFLOBDD.MkAnd(alpha1, beta)

    for i in range(dim//2, dim):
        for j in range(dim//2):
            if i == 0 and j == 0:
                continue
            else:
                nbrs = get_neighbors(i, j, dim)
                beta = pwc.CFLOBDD.MkFalse(level)
                for a in range(len(nbrs)):
                    a_idx = get_idx(nbrs[a][0], nbrs[a][1], dim)
                    gamma = pwc.CFLOBDD.MkTrue(level)
                    for b in range(len(nbrs)):
                        if a != b:
                            id_b = get_idx(nbrs[b][0], nbrs[b][1], dim)
                            gamma = pwc.CFLOBDD.MkAnd(gamma, pwc.CFLOBDD.MkNot(var_cflobdds[id_b-1]))

                    gamma = pwc.CFLOBDD.MkAnd(gamma, var_cflobdds[a_idx-1])
                    beta = pwc.CFLOBDD.MkOr(beta, gamma)

                idx = get_idx(i, j, dim)-1
                #beta = pwc.CFLOBDD.MkAnd(beta, var_cflobdds[idx])
                beta = pwc.CFLOBDD.MkOr(beta, pwc.CFLOBDD.MkNot(var_cflobdds[idx]))
                #list_of_betas.append(beta)
                alpha2 = pwc.CFLOBDD.MkAnd(alpha2, beta)
    
    
    for i in range(dim//2):
        for j in range(dim//2, dim):
            if i == 0 and j == 0:
                continue
            else:
                nbrs = get_neighbors(i, j, dim)
                beta = pwc.CFLOBDD.MkFalse(level)
                for a in range(len(nbrs)):
                    a_idx = get_idx(nbrs[a][0], nbrs[a][1], dim)
                    gamma = pwc.CFLOBDD.MkTrue(level)
                    for b in range(len(nbrs)):
                        if a != b:
                            id_b = get_idx(nbrs[b][0], nbrs[b][1], dim)
                            gamma = pwc.CFLOBDD.MkAnd(gamma, pwc.CFLOBDD.MkNot(var_cflobdds[id_b-1]))

                    gamma = pwc.CFLOBDD.MkAnd(gamma, var_cflobdds[a_idx-1])
                    beta = pwc.CFLOBDD.MkOr(beta, gamma)

                idx = get_idx(i, j, dim)-1
                #beta = pwc.CFLOBDD.MkAnd(beta, var_cflobdds[idx])
                beta = pwc.CFLOBDD.MkOr(beta, pwc.CFLOBDD.MkNot(var_cflobdds[idx]))
                #list_of_betas.append(beta)
                alpha3 = pwc.CFLOBDD.MkAnd(alpha3, beta)


    for i in range(dim//2,dim):
        for j in range(dim//2,dim):
            if i == 0 and j == 0:
                continue
            else:
                nbrs = get_neighbors(i, j, dim)
                beta = pwc.CFLOBDD.MkFalse(level)
                for a in range(len(nbrs)):
                    a_idx = get_idx(nbrs[a][0], nbrs[a][1], dim)
                    gamma = pwc.CFLOBDD.MkTrue(level)
                    for b in range(len(nbrs)):
                        if a != b:
                            id_b = get_idx(nbrs[b][0], nbrs[b][1], dim)
                            gamma = pwc.CFLOBDD.MkAnd(gamma, pwc.CFLOBDD.MkNot(var_cflobdds[id_b-1]))

                    gamma = pwc.CFLOBDD.MkAnd(gamma, var_cflobdds[a_idx-1])
                    beta = pwc.CFLOBDD.MkOr(beta, gamma)

                idx = get_idx(i, j, dim)-1
                #beta = pwc.CFLOBDD.MkAnd(beta, var_cflobdds[idx])
                beta = pwc.CFLOBDD.MkOr(beta, pwc.CFLOBDD.MkNot(var_cflobdds[idx]))
                #list_of_betas.append(beta)
                alpha4 = pwc.CFLOBDD.MkAnd(alpha4, beta)

    alpha = pwc.CFLOBDD.MkAnd(pwc.CFLOBDD.MkAnd(alpha1, alpha2),pwc.CFLOBDD.MkAnd(alpha3, alpha4))
    return alpha



pwc.CFLTests.init()
N = 2
alpha = build_constraints(N)
pwc.CFLOBDD.Print(alpha)
probs = [0.05 * (i + 1) for i in range(N)]
start = time.time()
sl = pwc.CFLOBDD.compute_prob(alpha, probs)
end = time.time()
print(sl)
print("time: ", (end - start))
pwc.CFLTests.clear()
