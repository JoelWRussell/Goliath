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
;;; The cross over function must take two polynomials and produce two child polynomials. 
;; **

;; @@
; Random number generator (excluding 0)

(defn random-num
  []
  "Generates a random number excluding the value 0"
  (let [x (rand)]
    (if (= x 0)
       (+ x (rand))) 
      x))


;need to add a final check if x is 0.
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/random-num</span>","value":"#'goliath/random-num"}
;; <=

;; @@
;Take to vectors    [[0 2]] [[2 3] [2 3]]]    [[0 2] [1 3] [3 5]]
; Created two empty chidren of length equal to partents.
; Assign genes from parent one to child one with probability X
; Else assign gene from parent one to child to poraility 1-X
; So long as the child genes do not exceed thier parental legnths add genes from parent 1 to offspring with the same probaiblity of.



(defn splice 
  [parent-gene child]
  (conj child parent-gene))


(defn cross-over [poly1 poly2]
  (let [c1 []
    	c2 []]  
(mapv (fn [poly1 poly2] (conj [(vec x)] (vec y))) poly1 poly2)
   ))

(contains? [[1 3] [1 2]] [1 3])


; the first genotype from parent 1 goes to child 1 with probability x
; else genotype one goes to child 2.
; offspring should have the same length as their parents.
; distribute polynimal 1 first



;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-unkown'>false</span>","value":"false"}
;; <=

;; @@
(cross-over [1 2 3] [1 2 3])
;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-long'>1</span>","value":"1"},{"type":"html","content":"<span class='clj-long'>2</span>","value":"2"},{"type":"html","content":"<span class='clj-long'>3</span>","value":"3"}],"value":"[1 2 3]"}],"value":"[[1 2 3]]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-long'>1</span>","value":"1"},{"type":"html","content":"<span class='clj-long'>2</span>","value":"2"},{"type":"html","content":"<span class='clj-long'>3</span>","value":"3"}],"value":"[1 2 3]"}],"value":"[[1 2 3]]"}],"value":"[[[1 2 3]] [[1 2 3]]]"}
;; <=

;; @@
(def x [[1 3] [3 1]])
(def y [[2 4] [1 0]])
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/y</span>","value":"#'goliath/y"}
;; <=

;; @@
(mapv (fn [x y] (conj [(vec x)] (vec y))) x y)

;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-long'>1</span>","value":"1"},{"type":"html","content":"<span class='clj-long'>3</span>","value":"3"}],"value":"[1 3]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-long'>2</span>","value":"2"},{"type":"html","content":"<span class='clj-long'>4</span>","value":"4"}],"value":"[2 4]"}],"value":"[[1 3] [2 4]]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-long'>3</span>","value":"3"},{"type":"html","content":"<span class='clj-long'>1</span>","value":"1"}],"value":"[3 1]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-long'>1</span>","value":"1"},{"type":"html","content":"<span class='clj-long'>0</span>","value":"0"}],"value":"[1 0]"}],"value":"[[3 1] [1 0]]"}],"value":"[[[1 3] [2 4]] [[3 1] [1 0]]]"}
;; <=

;; @@

;; @@
