
# Timsquery

## Where are we in the life cycle?

- The library is in a very early stage of development.
    - I cannot assure stability or bug-freeness.
    - We are still deciding what the API should be and what the scope of the project overall is.

1. Push to main.
2. Branched but fast moving <- **We are here**
3. Stable api but features might be dropped without notice.
4. Stable api and deprecations on changes.

## What is this?

Timsquery is meant to be a library and command line tool that allows querying TIMS data in a generic way.
Basically, pick a way in which you want your data to be aggregated, pick how you want to query it and pick
your file, and you get back results that match those three things!

More explicitly:
- The main design is to have modular components:
    - aggregators
    - indices
    - queries
    - tolerances

Thus, depending on the access pattern and purpose, the data can be queried in
different ways and aggregated in different ways (if you need random
access, use the index that works for that, if you need bulk
sequential, use that).

## What does the cli look like right now?


##  TODO:

- Add logging levels to instrumentations.
- Add missing_docks_in_private_items to clippy.
- Implement predicate pushdown on the setup of indices.
- Implement predicate pushdown on query execution for raw index.
