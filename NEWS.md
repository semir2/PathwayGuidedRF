# PathwayGuidedRF 0.5.0
major changes:
- changed mtry parameter to mtry.prop and mtry.pow to either apply mtry as power 
or proportion.
- included a downsampling procedure which generates the inbag samples for each 
tree by balanced resampling using the size of the minority class.


# PathwayGuidedRF 0.4.0

major changes:
- use more efficient permutation test in pw.rf.pred.error
  (which replaces adaptive test which includes a bug)
- added argument p.adjust.method to pw.rf.hunt, pw.rf.lefe and pw.rf.pred.error
  to control type of multiple testing adjustment of P values. (Note that the 
  default changed from the Bonferroni to the Benjamini-Hochberg procedure).

# PathwayGuidedRF 0.3.1

minor changes:
- allow underscore in variable names
- allow to specify number of cases and controls in function sim.data.study.2
- return gene specific importance in function pw.rf.hunt
  
bug fixes:
- enable one element vector of effects in function sim.data.study.2

# PathwayGuidedRF 0.3.0

initial commit to GitHub
