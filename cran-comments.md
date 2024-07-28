## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

* checking examples ... [32s] NOTE
  Examples with CPU (user + system) or elapsed time > 5s
              user system elapsed
  capzone    16.41   0.28   17.02
  tracelines 12.42   0.33   14.34
  
  The underlying computations require numerical integration which is slow. 
  Reducing the run time of the examples would require changing them in such a way that 
  they are no longer meaningful examples.
  
