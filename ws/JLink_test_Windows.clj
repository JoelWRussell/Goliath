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

;; @@
(def kernel (com.wolfram.jlink.MathLinkFactory/createKernelLink "-linkmode launch -linkname \"c:/program files/wolfram research/mathematica/10.0/mathkernel.exe\""))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;warm-swamp/kernel</span>","value":"#'warm-swamp/kernel"}
;; <=

;; @@
(.waitForAnswer kernel)
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-unkown'>8</span>","value":"8"}
;; <=

;; @@
(.getString kernel)
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-string'>&quot;In[1]:= &quot;</span>","value":"\"In[1]:= \""}
;; <=

;; @@
(.evaluate kernel "1 + 1")
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; @@
(.waitForAnswer kernel)
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-unkown'>3</span>","value":"3"}
;; <=

;; @@
(.getString kernel)
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-string'>&quot;2&quot;</span>","value":"\"2\""}
;; <=

;; @@
(.close kernel)
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; **
;;; This tests the score function directly. Note that I've changed the `InitFunctions` method so that it takes: the path to the MathKernel executable; the path to the directory that contains the package and the data; and the name of the data file. This makes it a bit more flexible.
;; **

;; @@
(goliath.mathlink.LagrangianScore/InitFunctions
  "c:/program files/wolfram research/mathematica/10.0/mathkernel.exe"
  "resources/"
  "mma_double.csv")
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; @@
(into [] (goliath.mathlink.LagrangianScore/GetScore (into-array Integer/TYPE [1, 1, 1, 1])))
;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>14.072445026423615</span>","value":"14.072445026423615"},{"type":"html","content":"<span class='clj-double'>8.129417044365267</span>","value":"8.129417044365267"}],"value":"[14.072445026423615 8.129417044365267]"}
;; <=

;; **
;;; Benchmark the performance of the score function with criterium.
;; **

;; @@
(defn test-score
  [ind]
  (into [] (goliath.mathlink.LagrangianScore/GetScore (into-array Integer/TYPE (flatten ind)))))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;warm-swamp/test-score</span>","value":"#'warm-swamp/test-score"}
;; <=

;; @@
(criterium/quick-bench (test-score [1, 1, 1, 1]))
;; @@
;; ->
;;; WARNING: Final GC required 3.462162560564387 % of runtime
;;; WARNING: Final GC required 28.43879097609376 % of runtime
;;; Evaluation count : 42 in 6 samples of 7 calls.
;;;              Execution time mean : 16.109125 ms
;;;     Execution time std-deviation : 720.667129 µs
;;;    Execution time lower quantile : 15.594178 ms ( 2.5%)
;;;    Execution time upper quantile : 17.336672 ms (97.5%)
;;;                    Overhead used : 2.718398 ns
;;; 
;;; Found 1 outliers in 6 samples (16.6667 %)
;;; 	low-severe	 1 (16.6667 %)
;;;  Variance from outliers : 13.8889 % Variance is moderately inflated by outliers
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
