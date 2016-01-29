;; gorilla-repl.fileformat = 1

;; **
;;; #Golliath 
;; **

;; @@
(ns goliath
  (:require [gorilla-plot.core :as plot]
            [regression-test-score :as ts]
            [clojure.repl :as repl]
            [clojure.pprint :as pprint]
            [clojure.walk :as walk]
            [clojure.zip :as zip]
            [gorilla-repl.table :as table]
            [gorilla-repl.html :as html]
    		[darwin.core :as darwin]
            [darwin.evolution.metrics :as metrics]
            [algebolic.expression.core :as expression]
            [algebolic.expression.tree :as tree]
            [algebolic.expression.genetics :as genetics]
            [algebolic.expression.score :as score]
            [algebolic.expression.render :as spea-render]
            [algebolic.expression.interpreter :as interpreter]
            [darwin.evolution.core :as evolution]
            [darwin.evolution.metrics :as metrics]
            [darwin.evolution.reproduction :as reproduction]
            [darwin.evolution.scoring :as scoring]
            [darwin.evolution.selection :as selection]
            [darwin.evolution.transform :as transform]
            [darwin.evolution.pareto :as pareto]
            [darwin.algorithms.spea2 :as spea2]
            [criterium.core :as criterium]))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; **
;;; **Genetic-algorithm declared variables required in the run**
;; **

;; @@
(def prob_inheritance 0.75)   ;probability that terms are directly inherited from parent to child.
;; @@

;; **
;;; ##Create New Individual 
;;; 
;; **

;; **
;;; 
;;; ##Crossover function
;;; 
;; **

;; **
;;; Helper function to distribute genes of a parent according to a probability r ∈ {0,1}
;; **

;; @@
(defn distribute-parent   
  "Distributes the alleses of a parent to children with seeded probability"
 ([parent-1 prob]
  (let [{c1 true c2 false} (group-by (fn [x] (<= prob (rand))) parent-1)]    
[c1 c2]  
   ))
  )
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/distribute-parent</span>","value":"#'goliath/distribute-parent"}
;; <=

;; **
;;; Cross-over function
;; **

;; @@
(defn cross-over
  "takes two individuals and crosses over the polynomial from one parent to another under a weighted selection bias" 
  [parent-1 parent-2] 
  (let [offspring (mapv vec (mapv set (mapv concat (distribute-parent (:genotype parent-1) prob_inheritance) (reverse (distribute-parent (:genotype parent-2) prob_inheritance)))))]
    [(offspring 0) (offspring 1)]
    ))
;; @@

;; @@
(cross-over (create-genotype 6 4 2) (create-genotype 6 4 2) 0.75)
;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>2</span>","value":"2"},{"type":"html","content":"<span class='clj-unkown'>4</span>","value":"4"}],"value":"[2 4]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>4</span>","value":"4"},{"type":"html","content":"<span class='clj-unkown'>4</span>","value":"4"}],"value":"[4 4]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>0</span>","value":"0"},{"type":"html","content":"<span class='clj-unkown'>1</span>","value":"1"}],"value":"[0 1]"}],"value":"[[2 4] [4 4] [0 1]]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>3</span>","value":"3"},{"type":"html","content":"<span class='clj-unkown'>0</span>","value":"0"}],"value":"[3 0]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>1</span>","value":"1"},{"type":"html","content":"<span class='clj-unkown'>3</span>","value":"3"}],"value":"[1 3]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>3</span>","value":"3"},{"type":"html","content":"<span class='clj-unkown'>1</span>","value":"1"}],"value":"[3 1]"}],"value":"[[3 0] [1 3] [3 1]]"}],"value":"[[[2 4] [4 4] [0 1]] [[3 0] [1 3] [3 1]]]"}
;; <=
