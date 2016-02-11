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
(def kernel (com.wolfram.jlink.MathLinkFactory/createKernelLink "-linkmode launch -linkname \"/Applications/Mathematica.app/Contents/MacOS/MathKernel\""))
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
  "/Applications/Mathematica.app/Contents/MacOS/MathKernel"
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
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>18.853345328976403</span>","value":"18.853345328976403"},{"type":"html","content":"<span class='clj-double'>-0.9999999995150077</span>","value":"-0.9999999995150077"}],"value":"[18.853345328976403 -0.9999999995150077]"}
;; <=

;; **
;;; Benchmark the performance of the score function with criterium.
;; **

;; @@
(defn test-score
  []
  (into [] (goliath.mathlink.LagrangianScore/GetScore (into-array Integer/TYPE [1, 1, 1, 1]))))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;warm-swamp/test-score</span>","value":"#'warm-swamp/test-score"}
;; <=

;; @@
(criterium/quick-bench (test-score))
;; @@
;; ->
;;; WARNING: Final GC required 21.2526281226028 % of runtime
;;; Evaluation count : 36 in 6 samples of 6 calls.
;;;              Execution time mean : 24.776276 ms
;;;     Execution time std-deviation : 5.432053 ms
;;;    Execution time lower quantile : 19.958330 ms ( 2.5%)
;;;    Execution time upper quantile : 31.245726 ms (97.5%)
;;;                    Overhead used : 1.882679 ns
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
