package org.trishullab.cflobdd.api.jdd;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import org.trishullab.cflobdd.CFLOBDD;
import jdd.bdd.OptimizedCache;
import jdd.bdd.SimpleCache;

class Configuration {
   public static int bddOpcacheDiv = 1;
   public static int bddNegcacheDiv = 2;
   public static int bddQuantcacheDiv = 3;
}

public class BDD {
   protected static final int CACHE_AND = 0;
   protected static final int CACHE_OR = 1;
   protected static final int CACHE_XOR = 2;
   protected static final int CACHE_BIIMP = 3;
   protected static final int CACHE_IMP = 4;
   protected static final int CACHE_NAND = 5;
   protected static final int CACHE_NOR = 6;
   protected static final int CACHE_RESTRICT = 7;
   protected static final int CACHE_EXISTS = 0;
   protected static final int CACHE_FORALL = 1;

    private final int level;
    private int varCount;
    private List<CFLOBDD> bddtable;
    private List<Integer> hashList;
    private List<Integer> refCount;
    protected SimpleCache op_cache;
    protected SimpleCache not_cache;
    protected SimpleCache quant_cache;

    /// A map from node's hash to node's index in bddtable
    private HashMap<Integer, Integer> lookup;

    private List<Integer> freeList;
    
    public BDD(int level, int cache_size) {
        this.level = level;
        bddtable = new ArrayList<>();
        refCount = new ArrayList<>();
        freeList = new LinkedList<>();
        hashList = new ArrayList<>();
        lookup = new HashMap<>();
        bddtable.add(CFLOBDD.createFalse(level));
        bddtable.add(CFLOBDD.createTrue(level));
        refCount.add(1);
        refCount.add(1);
        int hash0 = bddtable.get(0).hashCode();
        int hash1 = bddtable.get(1).hashCode();
        hashList.add(hash0);
        hashList.add(hash1);
        lookup.put(hash0, 0);
        lookup.put(hash1, 1);

        // from jdd

        this.op_cache = new OptimizedCache("OP", cache_size / Configuration.bddOpcacheDiv, 3, 2);
        this.not_cache = new OptimizedCache("NOT", cache_size / Configuration.bddNegcacheDiv, 1, 1);
        this.quant_cache = new OptimizedCache("QUANT", cache_size / Configuration.bddQuantcacheDiv, 3, 2);

    }

    public final int ref(int bdd) {
        bddtable.get(bdd).clone();
        refCount.set(bdd, refCount.get(bdd) + 1);
        // remove from freelist only if ref count is 1 (0 before)
        // if (refCount.get(bdd) == 1 && freeList.contains(bdd)) {
        //     System.err.println("Should not ref bdd in that is alreadly freed");
        //     // throw new RuntimeException("Should not ref bdd in that is alreadly freed");
        //     freeList.remove((Integer) bdd);
        //     lookup.put(hashList.get(bdd), bdd);
        // }
        return bdd;
    }
    public final int deref(int bdd) {
        if (bdd == 0 || bdd == 1) {
            return bdd;
        }
        if (refCount.get(bdd) == 0) {
            System.err.println("Double free");
            return bdd;
        }
        bddtable.get(bdd).deref();
        refCount.set(bdd, refCount.get(bdd) - 1);
        // don't free to enable caching.
        if (refCount.get(bdd) == 0) {
        // don't free to enable caching.
        //     freeList.add(bdd);
            lookup.remove(hashList.get(bdd));
        //     // table and hashlist are kept and
        }
        return bdd;
    }

    private int putBdd(CFLOBDD val) {
        int val_hash = val.hashCode();
        if (lookup.containsKey(val_hash)) {
            int ans = lookup.get(val_hash);
            val.deref();
            bddtable.get(ans).clone();
            refCount.set(ans, refCount.get(ans) + 1);
            return ans;
        }
        if (freeList.size() > 0) {
            int bdd = freeList.remove(0);
            bddtable.set(bdd, val);
            refCount.set(bdd, 1);
            lookup.put(val_hash, bdd);
            hashList.set(bdd, val_hash);
            return bdd;
        }
        bddtable.add(val);
        refCount.add(1);
        hashList.add(val_hash);
        lookup.put(val_hash, bddtable.size() - 1);
        return bddtable.size() - 1;
    }

    public int createVar() {
        CFLOBDD bdd = CFLOBDD.createVar(varCount++, level);
        return putBdd(bdd);
    }

    public int not(int bdd) {
        if (bdd < 2) {
            return bdd ^ 1; // no need to ref
        }
        if (this.not_cache.lookup(bdd)) {
            if (refCount.get(this.not_cache.answer) > 0) {
                return ref(this.not_cache.answer);
            } // otherwise it is invalidated
        }
        int hash = this.not_cache.hash_value;
        CFLOBDD notBdd = CFLOBDD.not(bddtable.get(bdd));
        int ans = putBdd(notBdd);
        this.not_cache.insert(hash, bdd, ans);
        return ans;
    }
        
    public int and(int bdd1, int bdd2) {
        if (bdd1 == 0 || bdd2 == 0) {
            return 0;
        }
        if (bdd1 == 1) {
            return ref(bdd2);
        }
        if (bdd2 == 1) {
            return ref(bdd1);
        }
        if (this.op_cache.lookup(bdd1, bdd2, CACHE_AND)) {
            if (refCount.get(this.op_cache.answer) > 0) {
                return ref(this.op_cache.answer);
            }
        }
        int hash = this.op_cache.hash_value;
        CFLOBDD andBdd = CFLOBDD.and(bddtable.get(bdd1), bddtable.get(bdd2));
        int ans = putBdd(andBdd);
        this.op_cache.insert(hash, bdd1, bdd2, CACHE_AND, ans);
        return ans;
    }
        
    public int or(int bdd1, int bdd2) {
        if (bdd1 == 1 || bdd2 == 1) {
            return 1;
        }
        if (bdd1 == 0) {
            return ref(bdd2);
        }
        if (bdd2 == 0) {
            return ref(bdd1);
        }
        if (this.op_cache.lookup(bdd1, bdd2, CACHE_OR)) {
            if (refCount.get(this.op_cache.answer) > 0) {
                return ref(this.op_cache.answer);
            }
        }
        int hash = this.op_cache.hash_value;
        CFLOBDD orBdd = CFLOBDD.or(bddtable.get(bdd1), bddtable.get(bdd2));
        int ans = putBdd(orBdd);
        this.op_cache.insert(hash, bdd1, bdd2, CACHE_OR, ans);
        return ans;
    }

    public int imp(int bdd1, int bdd2) {
        if (bdd1 == 0 || bdd2 == 1) {
            return 1;
        }
        if (bdd1 == 1) {
            return ref(bdd2);
        }
        if (bdd2 == 0) {
            return not(bdd1);
        }
        if (this.op_cache.lookup(bdd1, bdd2, CACHE_IMP)) {
            if (refCount.get(this.op_cache.answer) > 0) {
                return ref(this.op_cache.answer);
            }
        }
        int hash = this.op_cache.hash_value;
        CFLOBDD impBdd = CFLOBDD.imp(bddtable.get(bdd1), bddtable.get(bdd2));
        int ans = putBdd(impBdd);
        this.op_cache.insert(hash, bdd1, bdd2, CACHE_IMP, ans);
        return ans;
    }
    public int andTo(int bdd1, int bdd2) {
        int ans = and(bdd1, bdd2);
        deref(bdd1);
        return ans;
    }
    public int orTo(int bdd1, int bdd2) {
        int ans = or(bdd1, bdd2);
        deref(bdd1);
        return ans;
    }
    public int exists(int bdd, int var) {
        //TODO: exist cache
        CFLOBDD existsBdd = CFLOBDD.exists(bddtable.get(bdd), var);
        return putBdd(existsBdd);
    }
    public int satCount(int bdd) {
        return bddtable.get(bdd).satCount();
    }
    public long getMemoryUsage() {
        //TODO: implement
        System.err.println("not implemented");
        return 0;
    }
    public void oneSat(int bdd, int[] arr) {
        int[] ans = bddtable.get(bdd).oneSat();
        if (arr == null) {
            arr = ans;
        } else {
            for (int i = 0; i < arr.length; i++) {
                arr[i] = ans[i];
            }
        }
    }
    public boolean isValid(int bdd) {
        return bddtable.get(bdd).isValid();
    }

    public void printDot(String name, int bdd) {
        System.out.println("printDot: " + name);
        bddtable.get(bdd).print();
    }
    public int nodeCount(int bdd) {
        return bddtable.get(bdd).nodeCount();
    }
    public void save(int bdd, String path) {
        //TODO: implement
        System.err.println("not implemented");
    }
    public void cleanup() {
        CFLOBDD.dispose_module();
    }

    public static void main(String[] args) {
        BDD bdd = new BDD(10, 100000);
        int a = bdd.createVar();
        System.out.printf("a: %d, not a: %d\n", a, bdd.bddtable.get(a).hashCode());


        int a1 = bdd.createVar();
        int b1 = bdd.createVar();
        int c1 = bdd.createVar();
        int d1 = bdd.ref(bdd.and(a1, b1));
        int e1 = bdd.ref(bdd.and(d1, c1));

        int d2 = bdd.ref(bdd.and(b1, c1));
        int e2 = bdd.ref(bdd.and(d2, a1));
        System.out.printf("e1: %d, e2: %d\n", e1, e2);
        System.out.printf("e1: %x, e2: %x\n", bdd.bddtable.get(e1).getRoot(), bdd.bddtable.get(e2).getRoot());

        assert e1 == e2;
        bdd.cleanup();
    }
}
