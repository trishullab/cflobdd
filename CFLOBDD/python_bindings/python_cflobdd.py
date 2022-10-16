import python_wrapper_cflobdd as pwc

if __name__ == "__main__":

    pwc.CFLTests.init()
    x0 = pwc.CFLOBDD.MkProjection(0,2)
    x1 = pwc.CFLOBDD.MkProjection(1,2)
    x2 = pwc.CFLOBDD.MkNot(x0)
    a = pwc.CFLOBDD.MkAnd(x0,x1)
    a = pwc.CFLOBDD.MkOr(a, x2)
    probs = [0.5, 0.6]
    b = pwc.CFLOBDD.compute_prob(a, probs)
    print(b)
    pwc.CFLTests.clear()
