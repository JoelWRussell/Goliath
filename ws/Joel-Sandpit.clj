;; gorilla-repl.fileformat = 1

;; **
;;; # Joel-Sandpit
;;; 
;;; A place for trying out the cross-over and mutate functions.
;;; 
;;; The score function can be called with ts/poly-score
;;; 
;; **

;; @@
(ns goliath
  (:require [gorilla-plot.core :as plot]
            [score.regression-test-score :as ts]))
;; @@

;; @@
(ts/poly-score [1 3])
;; @@

;; **
;;; #Create New Individual 
;; **

;; @@
(defn new-poly-term
  [p_max n]
  "create a new polynomial term with maximum power p_max and maximum term length n_max"
  (vec (repeatedly n #(rand-int (+ 1 p_max))))
  )
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/new-poly-term</span>","value":"#'user/new-poly-term"}
;; <=

;; @@
(new-poly-term 5 4)
;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>4</span>","value":"4"},{"type":"html","content":"<span class='clj-unkown'>1</span>","value":"1"},{"type":"html","content":"<span class='clj-unkown'>5</span>","value":"5"},{"type":"html","content":"<span class='clj-unkown'>5</span>","value":"5"}],"value":"[4 1 5 5]"}
;; <=

;; @@
(defn new-individual
  [l p m]
  "creates a new individual of length between 1 and l_max decided at random.
  terms within an individual are of legnth n.
  polynomial terms are constrained to maximum power p_max"
  (let [l_max l
        p_max p 
        n m]
    (vec (repeatedly (+ 1 (rand-int (- l_max 1))) #(new-poly-term p_max n)))
    ))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/new-individual</span>","value":"#'user/new-individual"}
;; <=

;; @@
(new-individual 6 4 2)
;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>1</span>","value":"1"},{"type":"html","content":"<span class='clj-unkown'>4</span>","value":"4"}],"value":"[1 4]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>4</span>","value":"4"},{"type":"html","content":"<span class='clj-unkown'>0</span>","value":"0"}],"value":"[4 0]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>2</span>","value":"2"},{"type":"html","content":"<span class='clj-unkown'>0</span>","value":"0"}],"value":"[2 0]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>1</span>","value":"1"},{"type":"html","content":"<span class='clj-unkown'>0</span>","value":"0"}],"value":"[1 0]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>2</span>","value":"2"},{"type":"html","content":"<span class='clj-unkown'>0</span>","value":"0"}],"value":"[2 0]"}],"value":"[[1 4] [4 0] [2 0] [1 0] [2 0]]"}
;; <=

;; **
;;; 
;;; #Crossover function
;;; 
;; **

;; **
;;; Helper function to distribute genes of a parent according to a probability r ∈ {0,1}
;; **

;; @@
(defn distribute-parent   
 ([parent-1 prob]
  "Distributes the alleses of a parent to children with seeded probability"
  (let [{c1 true c2 false} (group-by (fn [x] (<= prob (rand))) parent-1)]    
[c1 c2]  
   ))
  )
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/distribute-parent</span>","value":"#'user/distribute-parent"}
;; <=

;; **
;;; **Construction tests - "distribute-parents"**
;; **

;; **
;;; #---------------
;; **

;; **
;;; Cross-over function
;; **

;; @@
(defn cross-over
  [parent-1 parent-2 prob] 
  "takes two individuals and crosses over the polynomial from one parent to another under a weighted selection bias" 
  (let [offspring (mapv vec (mapv set (mapv concat (distribute-parent parent-1 prob) (reverse (distribute-parent parent-2 prob)))))]
    [(offspring 0) (offspring 1)]
    )
    )
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/cross-over</span>","value":"#'user/cross-over"}
;; <=

;; @@
(cross-over (new-individual 6 4 2) (new-individual 6 4 2) 0.75)
;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>0</span>","value":"0"},{"type":"html","content":"<span class='clj-unkown'>0</span>","value":"0"}],"value":"[0 0]"}],"value":"[[0 0]]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>3</span>","value":"3"},{"type":"html","content":"<span class='clj-unkown'>4</span>","value":"4"}],"value":"[3 4]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>0</span>","value":"0"},{"type":"html","content":"<span class='clj-unkown'>3</span>","value":"3"}],"value":"[0 3]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>0</span>","value":"0"},{"type":"html","content":"<span class='clj-unkown'>1</span>","value":"1"}],"value":"[0 1]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>4</span>","value":"4"},{"type":"html","content":"<span class='clj-unkown'>0</span>","value":"0"}],"value":"[4 0]"}],"value":"[[3 4] [0 3] [0 1] [4 0]]"}],"value":"[[[0 0]] [[3 4] [0 3] [0 1] [4 0]]]"}
;; <=

;; @@
; prvent cross over producing nil
; prevent cross over producting a vector length greater than l_max
;; @@

;; **
;;; #Data Structure
;; **

;; **
;;; looking for a data strucutre of the form:  *( {:score -3 :genotype ((m1n1p1q1) (m2n2p2q2)...) } *
;;; 
;;; ------------------------- **Inputs** ------------------------- **Outputs**
;;; 
;;; - Crossover --------- 2 individuals ----------------- 2 Individuals
;;; - New Individual ------- nil -------------------------- Individual
;; **

;; @@
(defn create-individual
  (let [score]
    
  )
;; @@

;; @@
(defn cross
  []
  )
;; @@

;; @@
(new-individual 10 2 4)


m (count (first (new-individual 10 2 4)))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-unkown'>4</span>","value":"4"}
;; <=

;; @@

;; @@
