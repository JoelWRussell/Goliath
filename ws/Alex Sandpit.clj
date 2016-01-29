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
(def candidate {:score -3 :genotype [[0 4 5 1] [2 3 1 0]]})
(ts/ind-score candidate)
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-double'>35.35431659796228</span>","value":"35.35431659796228"}
;; <=

;; **
;;; ##Mutate function
;;; must take just one argument: the individual, access the genotype and mutate it, and then return the FULL individual.
;; **

;; @@
(def geno (:genotype candidate))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath.src/geno</span>","value":"#'goliath.src/geno"}
;; <=

;; **
;;; some helper funtions: rd-gene, a random gene generator and gene-replace, a function that accepts returns a vector 
;; **

;; @@
;;(defn rd-gene
;;  [powermax]
;;  [(rand-int powermax) (rand-int powermax)])

;; (rd-gene 8)

;; @@

;; @@
(defn rd-gene 
  [powermax spmax length]
  "where powermax is the maximum power of any given variable, spmax is the max sum of the powers and length is the length of the gene to be generated. spmax is checked and the function re-called if the gene is invalid"
  (let [new-gene (vec (repeatedly length #(rand-int (+ 1 powermax)))) sp (apply + new-gene)]
    (if (< sp spmax) new-gene
      (rd-gene powermax spmax length))
  )
 )


(rd-gene 10 10 4)
;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>0</span>","value":"0"},{"type":"html","content":"<span class='clj-unkown'>6</span>","value":"6"},{"type":"html","content":"<span class='clj-unkown'>3</span>","value":"3"},{"type":"html","content":"<span class='clj-unkown'>0</span>","value":"0"}],"value":"[0 6 3 0]"}
;; <=

;; @@
(count (first (:genotype {:score 5 :genotype [[3 5 5 6]]})))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-unkown'>4</span>","value":"4"}
;; <=

;; @@
(defn gene-replace
  [indv powermax spmax]
  "Replaces a gene from an individual (with a unique gene). In a given gene (term in the polynomial) no variable can have a higher power than powermax and the total powers of all variables cannot exceed spmax (sum power max)."
(let [length (count (first indv)) new-gene (rd-gene powermax spmax length) new-indv (assoc indv (rand-int (count indv)) new-gene) sp (apply + new-gene)]
  
  (if (= (count (set new-indv)) (count new-indv))
    new-indv
    (gene-replace indv powermax spmax)
 
 )))   
 
 ;;(cond
 ;;   (and (= (count (set new-indv)) (count new-indv)) (< sp spmax)) new-indv
 ;;   :else (gene-replace indv powermax spmax))
  



 (defn gene-add
  [indiv powermax spmax]
   "Adds a new gene to a GA individual. The gene will be unique addition to the individual. In a given gene (term in the polynomial) no variable can have a higher power than powermax and the total powers of all variables cannot exceed spmax (sum power max)."
  (let [length (count (first indiv)) new-gene (rd-gene powermax spmax length) new-indv (vec (conj (set indiv) new-gene)) sp (apply + new-gene)]
   
      (if (= (count new-indv) (count indiv))
      (gene-add indiv powermax spmax)
     new-indv
      )))


  ;;   (cond 
  ;;  (or (= (count new-indv) (count indiv)) (> sp spmax)) (gene-add indiv powermax spmax)
  ;;  :else new-indv)
  ;;  ))
 
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath.src/gene-add</span>","value":"#'goliath.src/gene-add"}
;; <=

;; @@
(gene-add geno 10 10)
;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-long'>2</span>","value":"2"},{"type":"html","content":"<span class='clj-long'>3</span>","value":"3"},{"type":"html","content":"<span class='clj-long'>1</span>","value":"1"},{"type":"html","content":"<span class='clj-long'>0</span>","value":"0"}],"value":"[2 3 1 0]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-long'>0</span>","value":"0"},{"type":"html","content":"<span class='clj-long'>4</span>","value":"4"},{"type":"html","content":"<span class='clj-long'>5</span>","value":"5"},{"type":"html","content":"<span class='clj-long'>1</span>","value":"1"}],"value":"[0 4 5 1]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>4</span>","value":"4"},{"type":"html","content":"<span class='clj-unkown'>0</span>","value":"0"},{"type":"html","content":"<span class='clj-unkown'>1</span>","value":"1"},{"type":"html","content":"<span class='clj-unkown'>4</span>","value":"4"}],"value":"[4 0 1 4]"}],"value":"[[2 3 1 0] [0 4 5 1] [4 0 1 4]]"}
;; <=

;; @@
(gene-replace geno 10 10)
;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>3</span>","value":"3"},{"type":"html","content":"<span class='clj-unkown'>2</span>","value":"2"},{"type":"html","content":"<span class='clj-unkown'>0</span>","value":"0"},{"type":"html","content":"<span class='clj-unkown'>1</span>","value":"1"}],"value":"[3 2 0 1]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-long'>2</span>","value":"2"},{"type":"html","content":"<span class='clj-long'>3</span>","value":"3"},{"type":"html","content":"<span class='clj-long'>1</span>","value":"1"},{"type":"html","content":"<span class='clj-long'>0</span>","value":"0"}],"value":"[2 3 1 0]"}],"value":"[[3 2 0 1] [2 3 1 0]]"}
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
