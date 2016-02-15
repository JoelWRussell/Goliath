;; gorilla-repl.fileformat = 1

;; **
;;; # Testing JLink
;;; 
;;; You can use this worksheet to test whether JLink is set up correctly on your computer.
;; **

;; @@
(ns warm-swamp
    (:require [gorilla-plot.core :as plot]
              [criterium.core :as criterium])
    (:import [com.wolfram.jlink]
      [goliath.mathlink]))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; **
;;; First, test the low level kernel link. This will reveal problems with the connection to MMA that are independent from the Goliath code.
;; **

;; **
;;; This tests the score function directly. Note that I've changed the `InitFunctions` method so that it takes: the path to the MathKernel executable; the path to the directory that contains the package and the data; and the name of the data file. This makes it a bit more flexible.
;; **

;; @@
(def mathKernelSz "c:/program files/wolfram research/mathematica/10.0/mathkernel.exe")
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;warm-swamp/mathKernelSz</span>","value":"#'warm-swamp/mathKernelSz"}
;; <=

;; @@
(def experimentalDataSz "mma_double.csv")
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;warm-swamp/experimentalDataSz</span>","value":"#'warm-swamp/experimentalDataSz"}
;; <=

;; @@
(def dt 0.1)
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;warm-swamp/dt</span>","value":"#'warm-swamp/dt"}
;; <=

;; @@

(def df 2)
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;warm-swamp/df</span>","value":"#'warm-swamp/df"}
;; <=

;; @@
(goliath.mathlink.LagrangianScore/InitFunctions
  mathKernelSz
  "resources/"
  experimentalDataSz
  0.1
  2
  )
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; @@
(into [] (goliath.mathlink.LagrangianScore/GetScore (into-array Integer/TYPE [1, 1, 1, 1]), 2))
;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>14.072445026423615</span>","value":"14.072445026423615"},{"type":"html","content":"<span class='clj-double'>8.129417044365265</span>","value":"8.129417044365265"}],"value":"[14.072445026423615 8.129417044365265]"}
;; <=

;; **
;;; Benchmark the performance of the score function with criterium.
;; **

;; @@
(defn test-score
  [ind]
  (into [] (goliath.mathlink.LagrangianScore/GetScore (into-array Integer/TYPE (flatten ind)) 2)))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;warm-swamp/test-score</span>","value":"#'warm-swamp/test-score"}
;; <=

;; @@
(criterium/quick-bench (test-score [1, 1, 1, 1]))
;; @@
;; ->
;;; WARNING: Final GC required 3.640581599655863 % of runtime
;;; WARNING: Final GC required 32.87504099568942 % of runtime
;;; Evaluation count : 42 in 6 samples of 7 calls.
;;;              Execution time mean : 15.310538 ms
;;;     Execution time std-deviation : 44.073140 µs
;;;    Execution time lower quantile : 15.251232 ms ( 2.5%)
;;;    Execution time upper quantile : 15.357332 ms (97.5%)
;;;                    Overhead used : 2.686370 ns
;;; 
;; <-
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; **
;;; Pretty good, 24ms per score function evaluation. Although it should be noted that there is some caching in the MMA kernel, so this might be artificially lowered.
;; **

;; @@
(goliath.mathlink.LagrangianScore/Shutdown)
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; @@

;; @@
