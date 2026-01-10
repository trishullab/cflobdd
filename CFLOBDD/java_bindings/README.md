# CFLOBDD Java Bindings

> [!WARNING]
> The Java binding is still under development; only part of the APIs are available.

## Introduction

This binding is a drop-in replacement of CFLOBDD for JDD. It uses `java.lang.foreign` to interface with the shared library of CFLOBDD, so Java 22 or higher is required.

The binding currently contains two types of APIs:
- `CFLOBDD.java`: the main class that provides the CFLOBDD APIs. It directly interacts with CFLOBDD shared library.
- `api/jdd/BDD.java`: is a wrapper around `CFLOBDD.java` that provides the same APIs as `jdd.bdd.BDD`. As it contains optimizations like operation cache, it is more efficient than `CFLOBDD.java`. **We recommend using `api/jdd/BDD.java` instead of `CFLOBDD.java` in most cases.**

## Usage

1. Replace `import jdd.bdd.BDD;` with `import org.trishullab.cflobdd.api.jdd.BDD;`
2. The constructor of `org.trishullab.cflobdd.api.jdd.BDD` is different from that of `jdd.bdd.BDD`. It takes `int level` instead of `int nodesize` as the first parameter. $n$-level CFLOBDD contains $2^n$ variables.
3. Make sure that `libcflobdd.so` or `libcflobdd.dylib` is in the `java.library.path`.