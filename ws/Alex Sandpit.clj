;; gorilla-repl.fileformat = 1

;; **
;;; # Alex Sandpit
;;; 
;;; A place for trying out the cross-over and mutate functions.
;;; 
;;; The individual-score-function can be called with ts/ind-score
;;; 
;; **

;; @@
(ns goliath.src
  (:require [gorilla-plot.core :as plot]
            [regression-test-score :as ts]))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; @@
(def candidate {:score -3 :genotype [[0 1] [2 0]]})
(ts/ind-score candidate)
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-double'>1.2112701606538574E-27</span>","value":"1.2112701606538574E-27"}
;; <=

;; **
;;; ##Mutate function
;;; must take just one argument: the individual, access the genotype and mutate it, and then return the FULL individual.
;; **

;; @@
(:genotype individual)
;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-long'>0</span>","value":"0"},{"type":"html","content":"<span class='clj-long'>1</span>","value":"1"}],"value":"[0 1]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-long'>2</span>","value":"2"},{"type":"html","content":"<span class='clj-long'>0</span>","value":"0"}],"value":"[2 0]"}],"value":"[[0 1] [2 0]]"}
;; <=

;; **
;;; some helper funtions: rd-gene, a random gene generator and gene-replace, a function that accepts returns a vector 
;; **

;; @@
(defn rd-gene
  [[ymax xmax]]
  [(rand-int ymax) (rand-int xmax)])

(rd-gene [3 5])

;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>2</span>","value":"2"},{"type":"html","content":"<span class='clj-unkown'>1</span>","value":"1"}],"value":"[2 1]"}
;; <=

;; @@
(defn gene-replace
  [indv [ymax xmax]]
(let [new-indv (assoc indv (rand-int (count indv)) (rd-gene [ymax xmax]))]
  (if (= (count (set new-indv)) (count new-indv))
    new-indv
    (gene-replace indv [ymax xmax])
    ))   
)
      
      
;;(defn gene-add
;;  [indv new-gene]
  
;;  (vec (conj (set indv) new-gene)))
 
(defn gene-add
  [individuall [ymax xmax]]
  (let [new-indv (vec (conj (set (:genotype individuall)) (rd-gene [ymax xmax])))
        indv (:genotype individuall)]
    (if (= (count new-indv) (count indv))
      (gene-add indv [ymax xmax])
     new-indv
      ))
  
  )
 
 
 
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath.src/gene-add</span>","value":"#'goliath.src/gene-add"}
;; <=

;; @@
(gene-add candidate [2 2] )
;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>1</span>","value":"1"},{"type":"html","content":"<span class='clj-unkown'>0</span>","value":"0"}],"value":"[1 0]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-long'>2</span>","value":"2"},{"type":"html","content":"<span class='clj-long'>0</span>","value":"0"}],"value":"[2 0]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-long'>0</span>","value":"0"},{"type":"html","content":"<span class='clj-long'>1</span>","value":"1"}],"value":"[0 1]"}],"value":"[[1 0] [2 0] [0 1]]"}
;; <=

;; @@
(gene-replace (:genotype individual) (rd-gene [3 3]))
;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>0</span>","value":"0"},{"type":"html","content":"<span class='clj-unkown'>1</span>","value":"1"}],"value":"[0 1]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-long'>2</span>","value":"2"},{"type":"html","content":"<span class='clj-long'>0</span>","value":"0"}],"value":"[2 0]"}],"value":"[[0 1] [2 0]]"}
;; <=

;; @@
(defn mutate
  [indv]
  (let [(:genotype indv) gen
        (
          k
          
          ) gen2]
        
       ;; {:score  }    
  
    
   
;; @@

;; @@

;; @@

;; @@

;; @@
