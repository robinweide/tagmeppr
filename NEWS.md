# tagmeppr 0.2

- [x] results-method to make a df of results
- [x] updated results-metadata to be more clear (orientation is now strand info)
- [x] `plotInsertions` now handles failed experiments without insertions
- [x] counts are now correct: `align()` did show some duplicate reads (which were not remove due to strand-info)
- [x] added the option to use padding in `findInsertions()`, which fixes issues with the mapper reporting matches over the TIS.
- [x] wrote checkPrimer to... check... the... primers...
- [x] `findInsertions()` now uses D_scores to set orientation as strand
- [x] `findInsertions()` uses the 5rev_3fwd-flag for setting the orientation
- [x] `plotInsertions()` has better spacing, colours and handles orientation.
- [x] `results()` now enables filtering on pvalue, counts and orientation
- [x] `plotSite()` is ready to use and also shows the orientation
- [x] a bug `plotSite()` regarding x-ticks is fixed
- [x] Write readme
- [x] checks if bwa and samtools are found
- [x] Removed the weird NA-warnings for emperical center in `align()`
- [x] PCR-duplicates are now removed on demand with `align(dedup = T)`
- [x] overhauled `findInsertions()`: handles failed experiments and outputs correct orientation

# tagmeppr 0.1

- [x] re-wrote indexing part to use internal ITRs
- [x] made classes for the index and the samples
- [x] packages is fully build-passed (without warnings) with travisCI
- [x] no longer need a fai: it is made internally w/ seqinfo
- [x] empirical center auto on 1K-Npad
- [x] correctly export methods
- [x] put postalign inside align
- [x] made a minimised fastq-seq to speed things up
- [x] made nice print methods
