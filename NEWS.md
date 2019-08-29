# tagmeppr 0.2

* results-method to make a df of results
* updated results-metadata to be more clear (orientation is now strand info)
* `plotInsertions` now handles failed experiments without insertions
* counts are now correct: `align()` did show some duplicate reads (which were not remove due to strand-info)
* added the option to use padding in `findInsertions()`, which fixes issues with the mapper reporting matches over the TIS.
* wrote checkPrimer to... check... the... primers...
* `findInsertions()` now uses D_scores to set orientation as strand
* `findInsertions()` uses the 5rev_3fwd-flag for setting the orientation
* `plotInsertions()` has better spacing, colours and handles orientation.
* `results()` now enables filtering on pvalue, counts and orientation
* `plotSite()` is ready to use and also shows the orientation
* a bug `plotSite()` regarding x-ticks is fixed

# tagmeppr 0.1

* re-wrote indexing part to use internal ITRs
* made classes for the index and the samples
* packages is fully build-passed (without warnings) with travisCI
* no longer need a fai: it is made internally w/ seqinfo
* empirical center auto on 1K-Npad
* correctly export methods
* put postalign inside align
* made a minimised fastq-seq to speed things up
* made nice print methods
